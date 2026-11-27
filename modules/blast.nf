#!/usr/bin/env nextflow

/*
 * BLAST any number of FASTA files in one pooled search: [sample, label, fasta] in, [sample, label, results] out.
 * Each search reads the whole database, so the time is set by the number of searches, not of reads: all files are
 * pooled, and a read is searched once per sample and cached (<store_dir>/blast_cache/<sample>).
 * Result lines are "all" or "nonhuman" (the search) then the BLAST columns, or "failed<TAB>read_id".
 */

include { blast_task; blast_search_pident } from './utils'

// The cache is outside the task directory, so containers need it mounted (created first, or docker creates it as root)
def blast_cache_root() {
    // one cache per search setting, so results from other settings are not reused
    def human = params.blast_human_db ? "_human-${file(params.blast_human_db).name}" : ""
    return file(params.store_dir).toAbsolutePath().normalize().toString() + "/blast_cache/${blast_task()}_pid${blast_search_pident()}${human}"
}

def blast_cache_mount() {
    if (!params.blast_cache || !workflow.containerEngine)
        return ''
    def dir = blast_cache_root()
    return workflow.containerEngine in ['singularity', 'apptainer'] ? "-B ${dir}" : "-v ${dir}:${dir}"
}

def blast_db_exists(dir, name) {
    return ["${name}.nal", "${name}.nsq", "${name}.00.nsq"].any{ file("${dir}/${it}").exists() }
}

// [directory, database name] for a database option: a path without extensions (/db/blast/core_nt), or a directory
// holding a database named after it (/db/blast/core_nt/core_nt.*) or one of the default names
def blast_db_parts(path = params.blast_db, default_names = ["core_nt", "nt"]) {
    def db = file(path)
    if (!db.isDirectory())
        return [db.parent, db.name]
    def name = ([db.name] + default_names).find{ blast_db_exists(db, it) }
    if (name)
        return [db, name]
    // a directory beside a database of the same name (/db/blast/core_nt/ and /db/blast/core_nt.*)
    return [db.parent, db.name]
}

// small for a remote search, so NCBI is not sent one huge request
def blast_chunk_reads() {
    return params.blast_chunk_reads ?: (params.blast_db ? 5000 : 20)
}

process pool_blast_queries {

    label "process_single"
    tag "all samples"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    containerOptions { blast_cache_mount() }
    publishDir path: { params.outdir }, mode: 'copy', pattern: "published/**", saveAs: { fn -> fn.replaceFirst(/^published\//, '') }

    input:
    // one "sample<TAB>label" line per file, in the order of the files
    val(manifest)
    path(queries, stageAs: "queries/query*")

    output:
    path("chunk_*.fasta"), emit: chunks, optional: true
    tuple path("map.tsv"), path("selected.tsv"), path("groups.tsv"), emit: tables
    path("published/*/*_blast_queries.fasta.gz"), emit: published

    script:
    cache = params.blast_cache ? blast_cache_root() : ""
    """
    cat > manifest.tsv <<'END_OF_MANIFEST'
${manifest.join('\n')}
END_OF_MANIFEST
    pool_blast_queries.py --manifest manifest.tsv --queries queries --cache "${cache}" --max_reads ${params.blast_max_reads ?: 0} --chunk_reads ${blast_chunk_reads()}
    """
}

// A remote search can fail or hang, so it is retried with a growing wait and timeout. A chunk that still fails is
// written as "failed" lines rather than failing the task, so the summaries are always produced.
process blastn_chunk {

    label "process_medium_plus_mem"
    label "process_long"
    tag "${chunk.baseName}"
    container "ncbi/blast:2.16.0"
    maxForks params.blast_db ? 4 : 1

    input:
    path(chunk)
    path(blast_db_dir)
    val(blast_db_name)
    path(human_db_dir, stageAs: "human_db_dir")
    val(human_db_name)

    output:
    path("results_${chunk.baseName}.txt")

    script:
    out_name = "results_${chunk.baseName}.txt"
    remote = !params.blast_db
    dbopts = remote ? "-db nt -remote" : "-db ${blast_db_dir}/${blast_db_name} -num_threads ${task.cpus}"
    outfmt = '"6 qseqid sacc sscinames staxids sstart send evalue pident length bitscore qcovs"'
    // the second search excludes human (9606) to find the best non-human hit, for leak tier 4
    limit_secs = remote ? (params.blast_remote_timeout as nextflow.util.Duration).toSeconds() : 24 * 3600
    search_opts = "-query ${chunk} -task ${blast_task()} -evalue 1e-6 -perc_identity ${blast_search_pident()} -max_target_seqs ${params.blast_max_targets} -outfmt ${outfmt}"
    blast_args = "${dbopts} ${search_opts}"
    search_all = "timeout \$((attempt * ${limit_secs})) blastn ${blast_args} -out all.tmp"
    search_nonhuman = "timeout \$((attempt * ${limit_secs})) blastn ${blast_args} -negative_taxids 9606 -out nonhuman.tmp"
    // with --blast_human_db, the best human hit comes from a human genome database (core_nt has no eukaryotic
    // chromosomes, so its human hits are clones and short records)
    human_db = params.blast_human_db ? true : false
    search_human = human_db ? "timeout \$((attempt * ${limit_secs})) blastn -db ${human_db_dir}/${human_db_name} -num_threads ${task.cpus} ${search_opts} -out human.tmp" : "true"
    result_files = (params.blast_nonhuman ? "all.txt nonhuman.txt" : "all.txt") + (human_db ? " human.txt" : "")
    wait_between = remote ? "sleep \$((attempt * 60))" : "true"
    """
    set -u
    ${remote ? "sleep 5" : "true"}
    ok=0
    for attempt in 1 2 3 4 5; do
      # a search which finished (renamed from .tmp) is not run again on a later attempt
      [ -e all.txt ] || { ${search_all} && mv all.tmp all.txt; } || true
      ${params.blast_nonhuman ? "[ -e nonhuman.txt ] || { ${search_nonhuman} && mv nonhuman.tmp nonhuman.txt; } || true" : "true"}
      ${human_db ? "[ -e human.txt ] || { ${search_human} && mv human.tmp human.txt; } || true" : "true"}
      if [ -e all.txt ] ${params.blast_nonhuman ? "&& [ -e nonhuman.txt ]" : ""} ${human_db ? "&& [ -e human.txt ]" : ""}; then ok=1; break; fi
      echo "BLAST attempt \$attempt failed" >&2
      ${wait_between}
    done
    : > result.tmp
    if [ "\$ok" = 1 ]; then
      for res in ${result_files}; do
        awk -v prefix="\${res%.txt}" '{ print prefix "\\t" \$0 }' "\$res" >> result.tmp
      done
    else
      grep '^>' ${chunk} | sed 's/^>//; s/[ \\t].*//' | awk '{ print "failed\\t" \$0 }' >> result.tmp || true
    fi
    mv result.tmp ${out_name}
    """
}

process split_blast_results {

    label "process_single"
    tag "all samples"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    containerOptions { blast_cache_mount() }
    publishDir path: { params.outdir }, mode: 'copy', saveAs: { fn -> fn.replaceFirst(/^out\//, '') }

    input:
    tuple path(map_file), path(selected), path(groups)
    path(results, stageAs: "results/*")

    output:
    path("out/*/blast/*.results_blastn.txt")

    script:
    cache = params.blast_cache ? blast_cache_root() : ""
    """
    mkdir -p results
    split_blast_results.py --map ${map_file} --selected ${selected} --groups ${groups} --results_dir results --cache "${cache}"
    """
}

workflow blast {
    take:
        query_ch    // [sample, label, fasta]

    main:
        if (params.blast_cache)
            file(blast_cache_root()).mkdirs()
        // a remote search needs no database, but the input must be a file
        db = params.blast_db ? blast_db_parts() : [file("${projectDir}/assets/NO_FILE", checkIfExists: true), "nt"]
        if (params.blast_human_db && !params.blast_db)
            error "--blast_human_db needs a local --blast_db for the other searches"
        human = params.blast_human_db ? blast_db_parts(params.blast_human_db, ["human_genome"]) : [file("${projectDir}/assets/NO_FILE", checkIfExists: true), "none"]
        if (params.blast_human_db && !blast_db_exists(human[0], human[1]))
            error "No BLAST database ${human[1]} in ${human[0]} (--blast_human_db ${params.blast_human_db}): give the database path without extensions"
        if (params.blast_db && !blast_db_exists(db[0], db[1]))
            error "No BLAST database ${db[1]} in ${db[0]} (--blast_db ${params.blast_db}): give the database path without extensions, e.g. /db/blast/core_nt for /db/blast/core_nt.*.nsq"

        query_ch
            .toList()
            .filter{ rows -> rows.size() > 0 }
            .multiMap{ rows ->
                manifest: rows.collect{ row -> "${row[0]}\t${row[1]}".toString() }
                files: rows.collect{ row -> row[2] }
            }
            .set{ pool_in }
        pool_blast_queries(pool_in.manifest, pool_in.files)
        blastn_chunk(pool_blast_queries.out.chunks.flatten(), db[0], db[1], human[0], human[1])
        split_blast_results(pool_blast_queries.out.tables, blastn_chunk.out.collect().ifEmpty([]))
        split_blast_results.out
            .flatten()
            .map{ f ->
                def sample = f.parent.parent.name
                [sample, f.name.substring(sample.length() + 1).replaceFirst(/\.results_blastn\.txt$/, ''), f]
            }
            .set{ results }

    emit:
        results
}

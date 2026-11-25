#!/usr/bin/env nextflow

/*
 * BLAST any number of FASTA files in one pooled search: [sample, label, fasta] in, [sample, label, results] out.
 * Each search reads the whole database, so the time is set by the number of searches, not of reads: all files are
 * pooled, and a read is searched once per sample and cached (<store_dir>/blast_cache/<sample>).
 * Result lines are "all" or "nonhuman" (the search) then the BLAST columns, or "failed<TAB>read_id".
 */

// The cache is outside the task directory, so containers need it mounted (created first, or docker creates it as root)
def blast_cache_root() {
    return file(params.store_dir).toAbsolutePath().normalize().toString() + "/blast_cache"
}

def blast_cache_mount() {
    if (!params.blast_cache || !workflow.containerEngine)
        return ''
    def dir = blast_cache_root()
    return workflow.containerEngine in ['singularity', 'apptainer'] ? "-B ${dir}" : "-v ${dir}:${dir}"
}

// [directory, database name] for --blast_db: a directory means the database nt in it
def blast_db_parts() {
    def db = file(params.blast_db)
    return db.isDirectory() ? [db, "nt"] : [db.parent, db.name]
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

    output:
    path("results_${chunk.baseName}.txt")

    script:
    out_name = "results_${chunk.baseName}.txt"
    remote = !params.blast_db
    dbopts = remote ? "-db nt -remote" : "-db ${blast_db_dir}/${blast_db_name} -num_threads ${task.cpus}"
    outfmt = '"6 qseqid sacc sscinames staxids sstart send evalue pident length bitscore qcovs"'
    // the second search excludes human (9606) to find the best non-human hit, for leak tier 4
    limit_secs = remote ? (params.blast_remote_timeout as nextflow.util.Duration).toSeconds() : 24 * 3600
    blast_args = "-query ${chunk} ${dbopts} -evalue 1e-6 -perc_identity 85 -max_target_seqs ${params.blast_max_targets} -outfmt ${outfmt}"
    search_all = "timeout \$((attempt * ${limit_secs})) blastn ${blast_args} -out all.tmp"
    search_nonhuman = "timeout \$((attempt * ${limit_secs})) blastn ${blast_args} -negative_taxids 9606 -out nonhuman.tmp"
    result_files = params.blast_nonhuman ? "all.txt nonhuman.txt" : "all.txt"
    wait_between = remote ? "sleep \$((attempt * 60))" : "true"
    """
    set -u
    ${remote ? "sleep 5" : "true"}
    ok=0
    for attempt in 1 2 3 4 5; do
      # a search which finished (renamed from .tmp) is not run again on a later attempt
      [ -e all.txt ] || { ${search_all} && mv all.tmp all.txt; } || true
      ${params.blast_nonhuman ? "[ -e nonhuman.txt ] || { ${search_nonhuman} && mv nonhuman.tmp nonhuman.txt; } || true" : "true"}
      if [ -e all.txt ] ${params.blast_nonhuman ? "&& [ -e nonhuman.txt ]" : ""}; then ok=1; break; fi
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
        if (params.blast_db && !file("${db[0]}/${db[1]}.nal").exists() && !file("${db[0]}/${db[1]}.nsq").exists() && !file("${db[0]}/${db[1]}.00.nsq").exists())
            error "No BLAST database ${db[1]} in ${db[0]} (--blast_db ${params.blast_db}): give the database path without extensions, e.g. /db/blast/core_nt"

        query_ch
            .toList()
            .filter{ rows -> rows.size() > 0 }
            .multiMap{ rows ->
                manifest: rows.collect{ row -> "${row[0]}\t${row[1]}".toString() }
                files: rows.collect{ row -> row[2] }
            }
            .set{ pool_in }
        pool_blast_queries(pool_in.manifest, pool_in.files)
        blastn_chunk(pool_blast_queries.out.chunks.flatten(), db[0], db[1])
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

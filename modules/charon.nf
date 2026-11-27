#!/usr/bin/env nextflow
include { resolve_path; minimap2_microbial; minimap2_host; extract_blast_queries; is_paired; reads_prefix; reads_glob; publish_reads; minimap2_unclassified; alignment_table } from '../modules/utils'

// Charon holds the index in memory plus ~7 GB (peak 14.2 GB with the 7.1 GB nanopore index) and is given 80% of the
// request as --max-memory, which must stay above the peak
def charon_memory(db, attempt) {
    return new nextflow.util.MemoryUnit((db.size() * 1.5) as long + 8.GB.toBytes()) * attempt
}

// A prebuilt Charon index from Zenodo, fetched once (already decompressed and md5-checked) into
// <store_dir>/charon_indexes/<record or "latest">/charon_<size>[_host_only].<tech>.idx
process fetch_charon_index {

    label "process_single"
    label "error_retry"
    container params.charon_container
    storeDir { "${params.store_dir}/charon_indexes/${params.charon_index_record ?: 'latest'}" }

    input:
    val(tech)
    val(size)
    val(host_only)

    output:
    path("charon_${size}${host_only ? '_host_only' : ''}.${tech}.idx")

    script:
    def record = params.charon_index_record ? "--record ${params.charon_index_record}" : ""
    """
    charon fetch --tech ${tech} --size ${size} ${host_only ? '--host-only' : ''} ${record} -o .
    """

    stub:
    """
    echo "${tech} ${size} ${host_only} ${params.charon_index_record}" > charon_${size}${host_only ? '_host_only' : ''}.${tech}.idx
    """
}

// The Charon index for a run: --db (or --db_host_only) if given, else the prebuilt one for --read_type from Zenodo.
// Host-only indexes are only built big.
workflow charon_index {
    take:
    host_only

    main:
    def db = host_only ? params.db_host_only : params.db
    if (db) {
        index = channel.value(file(resolve_path(db), type: "file", checkIfExists: true))
    } else {
        def tech = params.read_type == "illumina" ? "illumina" : "nanopore"
        index = fetch_charon_index(tech, host_only ? "big" : params.charon_index, host_only)
    }

    emit:
    index
}

// microbial_category is the category of reads kept: "microbial", or "non_host" for a host-only index
process run_charon {

    label "process_medium_plus_mem"
    label "process_long"
    memory { charon_memory(db, task.attempt) }
    tag "${unique_id}"
    container params.charon_container 
    maxForks 4

    publishDir path: { "${params.outdir}/${unique_id}/" }, mode: 'copy', pattern: "*.out"


    input:
    tuple val(unique_id), path(fastq)
    path(db)
    val(name)
    val(microbial_category)

    output:
    tuple val(unique_id), val(name), path("${reads_glob(reads_prefix(unique_id, fastq), microbial_category, fastq)}"), emit: microbial_fastq
    tuple val(unique_id), val(name), path("${reads_glob(reads_prefix(unique_id, fastq), 'human', fastq)}"), emit: human_fastq
    tuple val(unique_id), val(name), path("${unique_id}_${name}.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      ${params.charon_args} \
      --log charon_${unique_id}.log \
      --extract all \
      --prefix ${reads_prefix(unique_id, fastq)} \
      -t ${task.cpus} \
      --max-memory ${(task.memory.toGiga() * 0.8) as int} \
      > ${unique_id}_${name}.out

    """
}

// charon does not extract unclassified reads (status U in its report), so take them from the input by name
process extract_unclassified_reads {

    label "process_low"
    label "error_ignore"
    tag "${unique_id}"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"

    input:
    tuple val(unique_id), path(fastq), val(name), path(report)

    output:
    tuple val(unique_id), val(name), path("${unique_id}_${name}_unclassified*.fastq.gz")

    script:
    def files = fastq instanceof List ? fastq : [fastq]
    def extract = files.withIndex().collect{ f, i ->
        def out = "${unique_id}_${name}_unclassified${files.size() > 1 ? '_' + (i + 1) : ''}.fastq.gz"
        "zcat -f ${f} | awk 'NR==FNR {ids[\$1]; next} FNR%4==1 {id=substr(\$1, 2); sub(/\\/[12]\$/, \"\", id); keep=(id in ids)} keep' ids.txt - | gzip > ${out}"
    }
    """
    awk -F'\t' '\$1 == "U" {print \$2}' ${report} | sed 's#/[12]\$##' > ids.txt
    ${extract.join("\n    ")}
    """
}

// Run once per charon index; name labels it in all outputs
workflow evaluate_charon {
    take:
        fastq_ch
        db
        name
        microbial_category
        refs
    main:


    run_charon(fastq_ch, db, name, microbial_category)
    if ( params.publish_reads ) {
        extract_unclassified_reads(fastq_ch.combine(run_charon.out.result, by: 0))
        // one call for both sets of reads, since a process can only be called once
        publish_reads(run_charon.out.microbial_fastq.map{ id, method, reads -> [id, method, "microbial", reads] }
                        .mix(extract_unclassified_reads.out.map{ id, method, reads -> [id, method, "unclassified", reads] }))
        if ( params.map_unclassified ) {
            minimap2_unclassified(extract_unclassified_reads.out, refs)
            alignment_table(minimap2_unclassified.out)
        }
    }
    minimap2_microbial(run_charon.out.microbial_fastq, refs)
    minimap2_host(run_charon.out.human_fastq, refs)

    // BLAST runs once for all samples and methods, in evaluate_dehosting
    if ( params.evaluate_microbial ){
        extract_blast_queries(minimap2_microbial.out.chunk_sam_ch)
        query_ch = extract_blast_queries.out.queries
    } else {
        query_ch = channel.empty()
    }


    emit:
    report = run_charon.out.result
    microbial_sam = minimap2_microbial.out.sam_ch
    host_sam = minimap2_host.out
    blast_queries = query_ch
    no_reads = minimap2_microbial.out.no_reads


}

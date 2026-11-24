#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits; is_paired; reads_prefix; reads_glob; publish_reads; minimap2_unclassified; alignment_table } from '../modules/utils'

// Charon holds the whole index in memory and needs about 7 GB more than that for its reads and threads (peak_rss
// was 14.2 GB with the 7.1 GB nanopore index, and 4.3 GB with a 0.2 GB host-only index). It is given 80% of the
// request as --max-memory, which has to stay above the peak, so ask for 1.5 times the index plus 8 GB.
// A retry doubles it.
def charon_memory(db, attempt) {
    return new nextflow.util.MemoryUnit((db.size() * 1.5) as long + 8.GB.toBytes()) * attempt
}

// name labels the comparator (charon or charon_host_only) in all outputs. microbial_category is the extracted
// category holding the reads to keep: "microbial" for a normal index, "non_host" for a host-only index.
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
    // fastq is one file, or R1 and R2 of a pair (charon takes both and calls the pair together)
    """
    charon dehost ${fastq} \
      --db ${db} \
      --confidence 7 \
      --log charon_${unique_id}.log \
      --extract all \
      --prefix ${reads_prefix(unique_id, fastq)} \
      -t ${task.cpus} \
      --max-memory ${(task.memory.toGiga() * 0.8) as int} \
      > ${unique_id}_${name}.out

    """
}

process run_charon_chunk {

    label "process_medium_plus_mem"
    label "process_long"
    memory { charon_memory(db, task.attempt) }
    container params.charon_container
    maxForks 4

    publishDir path: { "${params.outdir}/${unique_id}/" }, mode: 'copy', pattern: "*.out"


    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("${fastq.baseName}_microbial.f*q.gz"), emit: microbial_fastq
    tuple val(unique_id), path("${fastq.baseName}_human.f*q.gz"), emit: human_fastq
    tuple val(unique_id), path("${unique_id}_charon.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      --confidence 7 \
      --log charon_${unique_id}.log \
      --extract all \
      --prefix ${fastq.baseName} \
      -t ${task.cpus} \
      --max-memory ${(task.memory.toGiga() * 0.8) as int} \
      > ${unique_id}_charon.out

    """
}

// charon does not extract the reads it leaves unclassified, but its report has a row (status U) for each of them, so
// take them from the input reads by name. For a pair, both reads of each pair.
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
    # status U is a read charon looked at but could not call
    awk -F'\t' '\$1 == "U" {print \$2}' ${report} | sed 's#/[12]\$##' > ids.txt
    ${extract.join("\n    ")}
    """
}

process cat_all_microbial_files {

    label "process_low"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"

    input:
    tuple val(unique_id), path(fastq_files)

    output:
    tuple val(unique_id), val("charon"), path("charon_${unique_id}_microbial.fq.gz")

    script:
    """
    cat \$(ls *.f*q*) > combined.tmp
    mv combined.tmp "charon_${unique_id}_microbial.fq.gz"
    """
}

process cat_all_human_files {

    label "process_low"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"

    input:
    tuple val(unique_id), path(fastq_files)

    output:
    tuple val(unique_id), val("charon"), path("charon_${unique_id}_human.fq.gz")

    script:
    """
    cat \$(ls *.f*q*) > combined.tmp
    mv combined.tmp "charon_${unique_id}_human.fq.gz"
    """
}

workflow run_charon_in_chunks {
    take:
    fastq_ch
    db

    main:
    fastq_ch.map{unique_id, fastq -> [unique_id, fastq.splitFastq(by: params.charon_chunk_size, file:true, compress:true)]}
            .transpose()
            .set{chunked_fastq_ch}

    run_charon_chunk(chunked_fastq_ch, db)
    
    run_charon_chunk.out.microbial_fastq.groupTuple() | cat_all_microbial_files | set{ microbial_fastq }

    run_charon_chunk.out.human_fastq.groupTuple() | cat_all_human_files | set{ human_fastq }

    run_charon_chunk.out.result.collectFile{ unique_id, result -> ["${unique_id}.charon.out", result.text]}
                                      .map { f -> [f.simpleName, "charon", f] }
                                      .set{ result }
    
    emit:
    microbial_fastq
    human_fastq
    result
}

// Run once per charon index: name labels the comparator in all outputs, and db is the index for it.
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
        // charon leaves some reads unclassified and does not extract them, so take them from the input by name
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

    if ( params.evaluate_microbial ){
        verify_microbial_host_hits(minimap2_microbial.out.chunk_sam_ch)
        // samples with no microbial reads have no BLAST step, so give them an empty result
        verify_microbial_host_hits.out
            .mix(minimap2_microbial.out.no_reads.map{ unique_id, method -> [unique_id, file("${projectDir}/resources/NO_BLAST")] })
            .set{ blast_ch }
    } else {
        blast_ch = Channel.empty()
    }


    emit:
    report = run_charon.out.result
    microbial_sam = minimap2_microbial.out.sam_ch
    host_sam = minimap2_host.out
    blast = blast_ch


}

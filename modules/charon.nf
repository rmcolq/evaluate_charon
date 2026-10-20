#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits } from '../modules/utils'

// name labels the comparator (charon or charon_host_only) in all outputs. microbial_category is the extracted
// category holding the reads to keep: "microbial" for a normal index, "non_host" for a host-only index.
process run_charon {

    label "process_medium_plus_mem"
    label "process_long"
    tag "${unique_id}"
    container { name == "charon_host_only" ? params.charon_host_only_container : params.charon_container }
    maxForks 4

    publishDir "${params.outdir}/${unique_id}/", mode: 'copy', pattern: "*.out"


    input:
    tuple val(unique_id), path(fastq)
    path(db)
    val(name)
    val(microbial_category)

    output:
    tuple val(unique_id), val(name), path("${fastq.baseName}_${microbial_category}.f*q.gz"), emit: microbial_fastq
    tuple val(unique_id), val(name), path("${fastq.baseName}_human.f*q.gz"), emit: human_fastq
    tuple val(unique_id), val(name), path("${unique_id}_${name}.out"),  emit: result

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
      > ${unique_id}_${name}.out

    """
}

process run_charon_chunk {

    label "process_medium_plus_mem"
    label "process_long"
    container 'docker.io/rmcolq/charon:v1.1.2'
    maxForks 4

    publishDir "${params.outdir}/${unique_id}/", mode: 'copy', pattern: "*.out"


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
    main:

    refs = file("$projectDir/${params.refs}", type: "file", checkIfExists:true)

    run_charon(fastq_ch, db, name, microbial_category)
    minimap2_microbial(run_charon.out.microbial_fastq, refs)
    minimap2_host(run_charon.out.human_fastq, refs)

    if ( params.evaluate_microbial ){
        verify_microbial_host_hits(minimap2_microbial.out.chunk_sam_ch)
        verify_microbial_host_hits.out.set{ blast_ch }
    } else {
        blast_ch = Channel.empty()
    }


    emit:
    report = run_charon.out.result
    microbial_sam = minimap2_microbial.out.sam_ch
    host_sam = minimap2_host.out
    blast = blast_ch


}

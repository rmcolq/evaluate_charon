#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits } from '../modules/utils'

process run_charon_chunk {

    label "process_medium_plus_mem"
    container 'docker.io/rmcolq/charon:v1.0.5'
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
      > ${unique_id}_charon.out

    """
}

process cat_all_microbial_files {

    label "process_low"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"

    input:
    tuple val(unique_id), path(fastq_files)

    output:
    tuple val(unique_id), val("charon"), path("charon_${unique_id}_microbial.f*q.gz")

    script:
    """
    cat \$(ls *.f*q*) >> "charon_${unique_id}_microbial.f*q.gz"
    """
}

process cat_all_human_files {

    label "process_low"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"

    input:
    tuple val(unique_id), path(fastq_files)

    output:
    tuple val(unique_id), val("charon"), path("charon_${unique_id}_human.f*q.gz")

    script:
    """
    cat \$(ls *.f*q*) >> "charon_${unique_id}_human.f*q.gz"
    """
}

workflow run_charon {
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

workflow evaluate_charon {
    take:
        fastq_ch
    main:

    db = file(params.db, type: "file", checkIfExists:true)
    refs = file("$projectDir/${params.refs}", type: "file", checkIfExists:true)

    run_charon(fastq_ch, db)
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

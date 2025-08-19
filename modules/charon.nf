#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits } from '../modules/utils'

process run_charon_chunk {

    label "process_medium"
    container 'docker.io/rmcolq/charon:v1.0.5'
    maxForks 4

    //publishDir "${params.outdir}/${unique_id}/intermediate/", mode: 'copy', pattern: "*.out"


    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), val("charon"), path("charon_${unique_id}_microbial.f*q"), emit: microbial_fastq
    tuple val(unique_id), val("charon"), path("charon_${unique_id}_human.f*q"), emit: human_fastq
    tuple val(unique_id), path("${unique_id}_charon.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      --confidence 7 \
      --log charon_${unique_id}.log \
      --extract all \
      --prefix charon_${unique_id} \
      -t ${task.cpus} \
      > ${unique_id}_charon.out

    PATTERN=(charon_${unique_id}*.f*q.gz)
    if [ -f \${PATTERN[0]} ]; then
        gunzip charon_${unique_id}*.f*q.gz
    fi
    """
}

process compress {
    input:
    tuple val(unique_id), val(method), path(fastq)

    output:
    tuple val(unique_id), val(method), path("${fastq}.gz")

    script:
    """
    gzip -c ${fastq} > "${fastq}.gz"
    """
}

process compress2 {
    input:
    tuple val(unique_id), val(method), path(fastq)

    output:
    tuple val(unique_id), val(method), path("${fastq}.gz")

    script:
    """
    gzip -c ${fastq} > "${fastq}.gz"
    """
}

workflow run_charon {
    take:
    fastq_ch
    db

    main:
    fastq_ch.map{unique_id, fastq -> [unique_id, fastq.splitFastq(by: params.chunk_size, file:true, compress:true)]}
            .transpose()
            .set{chunked_fastq_ch}
    run_charon_chunk(chunked_fastq_ch, db)
    
    run_charon_chunk.out.microbial_fastq.collectFile (storeDir: "${params.outdir}/intermediate"){ unique_id, method, result -> ["${unique_id}.charon_microbial.f*q", result.text] }
                                      .map { f -> [f.simpleName, "charon", f] }
                                      | compress
                                      | set{ microbial_fastq }

    run_charon_chunk.out.human_fastq.collectFile(storeDir: "${params.outdir}/intermediate"){ unique_id, method, result -> ["${unique_id}.charon_human.f*q", result.text] }
                                      .map { f -> [f.simpleName, "charon", f] }
                                      | compress2
                                      | set{ human_fastq }

    run_charon_chunk.out.result.collectFile(storeDir: "${params.outdir}/intermediate"){ unique_id, result -> ["${unique_id}.charon.out", result.text]}
                                      .map { f -> [f.simpleName, f] }
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

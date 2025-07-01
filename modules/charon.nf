#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; extract_microbial_host_hits; blastn_microbial_host_hits } from '../modules/utils'

process run_charon {

    label "process_medium_plus_mem"
    container 'docker.io/rmcolq/charon:v1.0.5'

    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("charon_${unique_id}_microbial.f*q.gz"), emit: microbial_fastq
    tuple val(unique_id), path("charon_${unique_id}_human.f*q.gz"), emit: human_fastq
    tuple val(unique_id), path("charon_${unique_id}.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      --confidence 7 \
      --log charon_${unique_id}.log \
      --extract all \
      --prefix charon_${unique_id} \
      -t ${task.cpus} \
      > charon_${unique_id}.out
    """
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
        blast_ch = Channel.empty()
        ref_bed = file("$projectDir/${params.ref_bed}", type: "file", checkIfExists:true)
        extract_microbial_host_hits(minimap2_microbial.out, ref_bed)
        if (params.blast_db)
            blast_db = file(params.blast_db, type: "dir", checkIfExists:true)
        else
            blast_db = file("$projectDir/${params.ref_bed}", type: "file", checkIfExists:true) // any file will do to not block
        blastn_microbial_host_hits(extract_microbial_host_hits.out, blast_db)
        blastn_microbial_host_hits.out.set{ blast_ch }
    } else {
        blast_ch = Channel.empty()
    }

    emit:
    report = run_charon.out.result
    microbial_sam = minimap2_microbial.out
    host_sam = minimap2_host.out
    blast = blast_ch


}

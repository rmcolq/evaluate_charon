#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits } from '../modules/utils'

process download_nohuman_index {
    label "process_single"
    storeDir "${params.store_dir}/nohuman/"
    container 'community.wave.seqera.io/library/nohuman:0.4.0--c0e60dcf9c349883'
    maxForks 1

    output:
        path("*.idx")

    script:
    """
    nohuman -d --db nohuman_v0.4.0.idx
    """
}
process run_nohuman {

    label "process_medium_plus_mem"
    label "process_long"
    container 'community.wave.seqera.io/library/nohuman:0.4.0--c0e60dcf9c349883'
    maxForks 2

    input:
    tuple val(unique_id), path(fastq)
    path(nohuman_index)

    output:
    tuple val(unique_id), val("nohuman"), path("nohuman_${unique_id}_microbial.fq.gz"), emit: microbial_fastq
    tuple val(unique_id), val("nohuman"), path("nohuman_${unique_id}_human.fq.gz"), emit: human_fastq
    tuple val(unique_id), path("nohuman_${unique_id}_microbial.fq.gz"), path("nohuman_${unique_id}_human.fq.gz"),  emit: combined

    script:
    """
    nohuman --db ${nohuman_index} -t ${task.cpus} -o "nohuman_${unique_id}_microbial.fq.gz" ${fastq}
    nohuman --human --db ${nohuman_index} -t ${task.cpus} -o "nohuman_${unique_id}_human.fq.gz" ${fastq}
    """
}

process collect_classifications {

    label "process_low"
    container 'community.wave.seqera.io/library/deacon:0.5.0--5e06f862e47bcb8a'

    publishDir "${params.outdir}/${unique_id}/", mode: 'copy', pattern: "*.out"

    input:
    tuple val(unique_id), path(microbial_fastq), path(host_fastq)

    output:
    tuple val(unique_id), val("nohuman"), path("${unique_id}.nohuman.out"),  emit: result

    script:
    """
    cat ${host_fastq}  | gunzip | awk 'NR%4==1 {print substr(\$1,2)}' > list_host
    cat ${microbial_fastq}  | gunzip | awk 'NR%4==1 {print substr(\$1,2)}' > list_microbial
    echo -e "read_id\tclassification" > "${unique_id}.nohuman.out"
    for id in \$(cat list_microbial)
      do
        echo -e "\$id\tmicrobial"
      done >> "${unique_id}.nohuman.out"
    for id in \$(cat list_host)
      do
        echo -e "\$id\thuman"
      done >> "${unique_id}.nohuman.out"
    """
}

workflow evaluate_nohuman {
    take:
        fastq_ch
    main:

    download_nohuman_index()
    refs = file("$projectDir/${params.refs}", type: "file", checkIfExists:true)

    run_nohuman(fastq_ch, download_nohuman_index.out)
    collect_classifications(run_nohuman.out.combined)
    minimap2_microbial(run_nohuman.out.microbial_fastq, refs)
    minimap2_host(run_nohuman.out.human_fastq, refs)

    if ( params.evaluate_microbial ){
        verify_microbial_host_hits(minimap2_microbial.out.chunk_sam_ch)
        verify_microbial_host_hits.out.set{ blast_ch }
    } else {
        blast_ch = Channel.empty()
    }

    emit:
        report = collect_classifications.out.result
        microbial_sam = minimap2_microbial.out.sam_ch
        host_sam = minimap2_host.out
        blast = blast_ch
}

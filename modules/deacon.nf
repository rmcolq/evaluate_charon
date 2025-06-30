#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; extract_microbial_host_hits; blastn_microbial_host_hits } from '../modules/utils'

process download_deacon_index {
    label "process_single"
    storeDir "${params.store_dir}/deacon/"
    output:
        path("*.idx")

    script:
    """
    wget ${params.deacon_index}
    """
}
process run_deacon {

    label "process_medium"
    container 'community.wave.seqera.io/library/deacon:0.5.0--5e06f862e47bcb8a'

    input:
    tuple val(unique_id), path(fastq)
    path(deacon_index)

    output:
    tuple val(unique_id), path("deacon_${unique_id}_microbial.fq.gz"), emit: microbial_fastq
    tuple val(unique_id), path("deacon_${unique_id}_human.fq.gz"), emit: human_fastq
    tuple val(unique_id), path("deacon_${unique_id}_microbial.fq.gz"), path("deacon_${unique_id}_human.fq.gz"),  emit: combined

    script:
    """
    deacon filter ${deacon_index} ${fastq} -o "deacon_${unique_id}_human.fq.gz"
    deacon filter -d ${deacon_index} ${fastq} -o "deacon_${unique_id}_microbial.fq.gz"
    """
}

process collect_classifications {

    label "process_low"

    input:
    tuple val(unique_id), path(microbial_fastq), path(host_fastq)

    output:
    tuple val(unique_id), path("deacon_${unique_id}.out"),  emit: result

    script:
    """
    zgrep ">" ${microbial_fastq} | cut -f1 -d" " | cut -f2 -d'>' > list_microbial
    zgrep ">" ${host_fastq} | cut -f1 -d" " | cut -f2 -d'>' > list_host
    #echo -e "status\tread_id\tclassification" > "deacon_${unique_id}.out"
    for id in \$(cat list_microbial)
      do
        echo -e "\$id\tmicrobial"
      done >> "deacon_${unique_id}.out"
    for id in \$(cat list_host)
      do
        echo -e "\$id\thuman"
      done >> "deacon_${unique_id}.out"
    """
}

workflow evaluate_deacon {
    take:
        fastq_ch
    main:

    download_deacon_index()
    refs = file("$projectDir/${params.refs}", type: "file", checkIfExists:true)

    run_deacon(fastq_ch, download_deacon_index.out)
    collect_classifications(run_deacon.out.combined)
    minimap2_microbial(run_deacon.out.microbial_fastq, refs)
    minimap2_host(run_deacon.out.human_fastq, refs)

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
        report = collect_classifications.out.result
        microbial_sam = minimap2_microbial.out
        host_sam = minimap2_host.out
        blast = blast_ch
}

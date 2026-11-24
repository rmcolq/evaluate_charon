#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits; publish_reads; is_paired; reads_prefix; reads_glob; first_file } from '../modules/utils'

process download_nohuman_index {
    label "process_single"
    storeDir "${params.store_dir}/nohuman/"
    container 'community.wave.seqera.io/library/nohuman:0.5.1--69df84d3602dde50'
    maxForks 1

    output:
        path("*.idx")

    script:
    """
    nohuman --download --db-version HPRC.r2 --db nohuman_v0.5.1.idx
    """
}
process run_nohuman {

    label "process_medium_plus_mem"
    tag "${unique_id}"
    label "process_long"
    container 'community.wave.seqera.io/library/nohuman:0.5.1--69df84d3602dde50'
    maxForks 2

    input:
    tuple val(unique_id), path(fastq)
    path(nohuman_index)

    output:
    tuple val(unique_id), val("nohuman"), path("${nohuman_out(unique_id, 'microbial', fastq)}"), emit: microbial_fastq
    tuple val(unique_id), val("nohuman"), path("${nohuman_out(unique_id, 'human', fastq)}"), emit: human_fastq
    tuple val(unique_id), path("${nohuman_out(unique_id, 'microbial', fastq)}"), path("${nohuman_out(unique_id, 'human', fastq)}"),  emit: combined

    script:
    // fastq is one file, or R1 and R2 of a pair
    """
    nohuman --db ${nohuman_index} -t ${task.cpus} ${nohuman_output_args(unique_id, 'microbial', fastq)} ${fastq}
    nohuman --human --db ${nohuman_index} -t ${task.cpus} ${nohuman_output_args(unique_id, 'human', fastq)} ${fastq}
    """
}

// nohuman options naming the output file(s) of a category: -o, or --out1/--out2 for a pair
def nohuman_output_args(unique_id, category, reads) {
    return is_paired(reads) ?
        "--out1 nohuman_${unique_id}_${category}_1.fq.gz --out2 nohuman_${unique_id}_${category}_2.fq.gz" :
        "-o nohuman_${unique_id}_${category}.fq.gz"
}

// Output file(s) of a category: one file, or _1 and _2 for a pair
def nohuman_out(unique_id, category, reads) {
    return is_paired(reads) ? "nohuman_${unique_id}_${category}_?.fq.gz" : "nohuman_${unique_id}_${category}.fq.gz"
}

process collect_classifications {

    label "process_low"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/nohuman:0.5.1--69df84d3602dde50'

    publishDir "${params.outdir}/${unique_id}/", mode: 'copy', pattern: "*.out"

    input:
    tuple val(unique_id), path(microbial_fastq), path(host_fastq)

    output:
    tuple val(unique_id), val("nohuman"), path("${unique_id}.nohuman.out"),  emit: result

    script:
    """
    # one id per read, or per pair (taken from R1)
    cat ${first_file(host_fastq)}  | gunzip | awk 'NR%4==1 {print substr(\$1,2)}' | sed 's#/[12]\$##' > list_host
    cat ${first_file(microbial_fastq)}  | gunzip | awk 'NR%4==1 {print substr(\$1,2)}' | sed 's#/[12]\$##' > list_microbial
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
        refs
    main:

    download_nohuman_index()

    run_nohuman(fastq_ch, download_nohuman_index.out)
    if ( params.publish_reads ) {
        publish_reads(run_nohuman.out.microbial_fastq.map{ id, method, reads -> [id, method, "microbial", reads] })
    }
    collect_classifications(run_nohuman.out.combined)
    minimap2_microbial(run_nohuman.out.microbial_fastq, refs)
    minimap2_host(run_nohuman.out.human_fastq, refs)

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
        report = collect_classifications.out.result
        microbial_sam = minimap2_microbial.out.sam_ch
        host_sam = minimap2_host.out
        blast = blast_ch
}

#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; verify_microbial_host_hits; mapping_reference; is_paired; reads_prefix; reads_glob; first_file } from '../modules/utils'

process download_deacon_index {
    label "process_single"
    storeDir "${params.store_dir}/deacon/"
    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'
    maxForks 1

    output:
        path("*.idx")

    script:
    """
    wget ${params.deacon_index}
    """
}
process run_deacon {

    label "process_medium"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'
    maxForks 2

    input:
    tuple val(unique_id), path(fastq)
    path(deacon_index)
    val(name)
    val(extra_args)

    output:
    tuple val(unique_id), val(name), path("${deacon_out(name, unique_id, 'microbial', fastq)}"), emit: microbial_fastq
    tuple val(unique_id), val(name), path("${deacon_out(name, unique_id, 'human', fastq)}"), emit: human_fastq
    tuple val(unique_id), val(name), path("${deacon_out(name, unique_id, 'microbial', fastq)}"), path("${deacon_out(name, unique_id, 'human', fastq)}"),  emit: combined

    script:
    // fastq is one file, or R1 and R2 of a pair (deacon counts a hit in either mate for the pair).
    // extra_args are applied to both calls so that the human and microbial sets are exact complements
    def outputs = { category ->
        is_paired(fastq) ?
            "-o ${name}_${unique_id}_${category}_1.fq.gz -O ${name}_${unique_id}_${category}_2.fq.gz" :
            "-o ${name}_${unique_id}_${category}.fq.gz"
    }
    """
    deacon filter ${extra_args} ${deacon_index} ${fastq} ${outputs('human')}
    deacon filter -d ${extra_args} ${deacon_index} ${fastq} ${outputs('microbial')}
    """
}

// Output file(s) of a category: one file, or _1 and _2 for a pair
def deacon_out(name, unique_id, category, reads) {
    return is_paired(reads) ? "${name}_${unique_id}_${category}_?.fq.gz" : "${name}_${unique_id}_${category}.fq.gz"
}

process collect_classifications {

    label "process_low"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/deacon:0.17.0--43cd5289edd1686c'

    publishDir "${params.outdir}/${unique_id}/", mode: 'copy', pattern: "*.out"

    input:
    tuple val(unique_id), val(name), path(microbial_fastq), path(host_fastq)

    output:
    tuple val(unique_id), val(name), path("${unique_id}.${name}.out"),  emit: result

    script:
    """
    # one id per read, or per pair (taken from R1)
    cat ${first_file(host_fastq)}  | gunzip | awk 'NR%4==1 {print substr(\$1,2)}' | sed 's#/[12]\$##' > list_host
    cat ${first_file(microbial_fastq)}  | gunzip | awk 'NR%4==1 {print substr(\$1,2)}' | sed 's#/[12]\$##' > list_microbial
    echo -e "read_id\tclassification" > "${unique_id}.${name}.out"
    for id in \$(cat list_microbial)
      do
        echo -e "\$id\tmicrobial"
      done >> "${unique_id}.${name}.out"
    for id in \$(cat list_host)
      do
        echo -e "\$id\thuman"
      done >> "${unique_id}.${name}.out"
    """
}

// Run once per deacon setting: name labels the comparator (e.g. deacon, deacon_sensitive) in all outputs
// and extra_args are added to the deacon filter calls.
workflow evaluate_deacon {
    take:
        fastq_ch
        deacon_index
        name
        extra_args
    main:

    refs = mapping_reference()

    run_deacon(fastq_ch, deacon_index, name, extra_args)
    collect_classifications(run_deacon.out.combined)
    minimap2_microbial(run_deacon.out.microbial_fastq, refs)
    minimap2_host(run_deacon.out.human_fastq, refs)

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

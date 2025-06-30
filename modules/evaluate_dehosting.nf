include { evaluate_charon } from '../modules/charon'
include { evaluate_deacon } from '../modules/deacon'

process evaluate_summary {

    label "process_low"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    publishDir "${params.outdir}/${unique_id}/", mode: 'copy'

    input:
    tuple val(unique_id), path(charon_report), path(host_sam), path(microbial_sam), path(blast_result), path(additional_report)

    output:
    path "${unique_id}_summary.csv", emit: summary
    path "${unique_id}_full.csv", emit: full
    path "${unique_id}*_data.csv", emit: data
    path "${unique_id}*_taxa.csv", emit: taxa, optional:true
    path "${unique_id}*_accs.csv", emit: accs, optional:true

    script:
    """
    evaluate.py \
      -i ${charon_report} \
      --microbial_sam ${microbial_sam} \
      --host_sam ${host_sam} \
      --blast_result ${blast_result} \
      --additional ${additional_report} \
      -p "charon_${unique_id}"
    """
}

workflow evaluate_dehosting {
    take:
    fastq_ch

    main:

    evaluate_charon(fastq_ch)
    evaluate_deacon(fastq_ch)

    evaluate_charon.out.host_sam.concat(evaluate_deacon.out.host_sam)
                                .collectFile( keepHeader:true, skip:29, newLine:true)
                                .set{ host_sam }

    evaluate_charon.out.microbial_sam.concat(evaluate_deacon.out.microbial_sam)
                                     .collectFile( keepHeader:true, skip:29, newLine:true)
                                     .set{ microbial_sam }

    evaluate_charon.out.blast.concat(evaluate_deacon.out.blast)
                             .collectFile(newLine:true)
                             .set{ blast }

    evaluate_charon.out.report
                 .combine(host_sam, by: 0)
                 .combine(microbial_sam, by: 0)
                 .combine(blast, by: 0)
                 .combine(evaluate_deacon.out.report, by:0)
                 .set{ eval_ch }

    evaluate_summary(eval_ch)
}




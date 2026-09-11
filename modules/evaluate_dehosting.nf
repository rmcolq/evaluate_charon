include { evaluate_charon } from '../modules/charon'
include { evaluate_deacon } from '../modules/deacon'
include { evaluate_nohuman } from '../modules/nohuman'
include { evaluate_summary } from '../modules/utils'




workflow evaluate_dehosting {
    take:
    fastq_ch

    main:

    if (params.charon){
      evaluate_charon(fastq_ch)
      evaluate_charon.out.report
                 .combine(evaluate_charon.out.host_sam, by: 0)
                 .combine(evaluate_charon.out.microbial_sam, by: 0)
                 .combine(evaluate_charon.out.blast, by: 0)
                 .view()
                 .set{ eval_ch }

      evaluate_summary(eval_ch)
    }

    if (params.deacon){
      evaluate_deacon(fastq_ch)
      evaluate_deacon.out.report
                 .combine(evaluate_deacon.out.host_sam, by: 0)
                 .combine(evaluate_deacon.out.microbial_sam, by: 0)
                 .combine(evaluate_deacon.out.blast, by: 0)
                 .view()
                 .set{ eval_ch }

      evaluate_summary(eval_ch)
    }

    if (params.nohuman){
      evaluate_nohuman(fastq_ch)
      evaluate_nohuman.out.report
                 .combine(evaluate_nohuman.out.host_sam, by: 0)
                 .combine(evaluate_nohuman.out.microbial_sam, by: 0)
                 .combine(evaluate_nohuman.out.blast, by: 0)
                 .view()
                 .set{ eval_ch }

      evaluate_summary(eval_ch)
    }

}




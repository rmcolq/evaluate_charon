include { evaluate_charon } from '../modules/charon'
include { evaluate_deacon } from '../modules/deacon'
include { evaluate_nohuman } from '../modules/nohuman'
include { evaluate_summary } from '../modules/utils'

// Combine the pieces for one method into the tuple evaluate_summary expects:
// [id, classifier, report, host_sam, microbial_sam, blast_result]
// When blast is not run, an empty placeholder is used (evaluate.py treats an empty file as no hits).
def build_eval_ch(method) {
    def no_blast = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def eval_ch = method.report
                    .combine(method.host_sam, by: 0)
                    .combine(method.microbial_sam, by: 0)
    if (params.evaluate_microbial) {
        eval_ch = eval_ch.combine(method.blast, by: 0)
    } else {
        eval_ch = eval_ch.map{ id, classifier, report, host, microbial -> [id, classifier, report, host, microbial, no_blast] }
    }
    return eval_ch
}

workflow evaluate_dehosting {
    take:
    fastq_ch

    main:
    eval_ch = Channel.empty()

    if (params.charon){
      evaluate_charon(fastq_ch)
      eval_ch = eval_ch.mix(build_eval_ch(evaluate_charon.out))
    }

    if (params.deacon){
      evaluate_deacon(fastq_ch)
      eval_ch = eval_ch.mix(build_eval_ch(evaluate_deacon.out))
    }

    if (params.nohuman){
      evaluate_nohuman(fastq_ch)
      eval_ch = eval_ch.mix(build_eval_ch(evaluate_nohuman.out))
    }

    // a process can only be invoked once per workflow, so all methods share this call
    contaminant_accs = params.contaminant_accs ? file(params.contaminant_accs, checkIfExists: true) : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    evaluate_summary(eval_ch, contaminant_accs)

    emit:
    done = evaluate_summary.out.done
}

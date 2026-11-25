include { evaluate_charon } from '../modules/charon'
include { evaluate_charon as evaluate_charon_host_only } from '../modules/charon'
include { evaluate_deacon; download_deacon_index } from '../modules/deacon'
include { evaluate_deacon as evaluate_deacon_sensitive } from '../modules/deacon'
include { evaluate_nohuman } from '../modules/nohuman'
include { evaluate_summary; mapping_index; resolve_path } from '../modules/utils'
include { blast } from '../modules/blast'

// [id, classifier, report, host_sam, host_sampled_ids, microbial_sam, blast_result] for evaluate_summary. A sample
// with no microbial reads, or a run without BLAST, gets an empty BLAST result.
def build_eval_ch(method, name, blast_results) {
    def no_blast = file("${projectDir}/resources/NO_BLAST", checkIfExists: true)
    def eval_ch = method.report
                    .combine(method.host_sam, by: 0)
                    .combine(method.microbial_sam, by: 0)
    if (params.evaluate_microbial) {
        def blast_ch = blast_results
                    .filter{ row -> row[1] == name }
                    .map{ id, m, result -> [id, result] }
                    .mix(method.no_reads.map{ id, m -> [id, no_blast] })
        eval_ch = eval_ch.combine(blast_ch, by: 0)
    } else {
        eval_ch = eval_ch.map{ id, classifier, report, host, host_ids, microbial -> [id, classifier, report, host, host_ids, microbial, no_blast] }
    }
    return eval_ch
}

// [id, classifier] once a method has classified and mapped a sample. Slots are freed on this, not on the summary,
// which waits for the pooled BLAST of all samples and so for every sample to pass the gate.
def reads_done(method) {
    return method.report
                .combine(method.host_sam, by: 0)
                .combine(method.microbial_sam, by: 0)
                .map{ row -> [row[0], row[1]] }
}

workflow evaluate_dehosting {
    take:
    fastq_ch

    main:
    query_ch = channel.empty()
    done_ch = channel.empty()

    mapping_index()
    refs = mapping_index.out.index

    if (params.charon){
      evaluate_charon(fastq_ch, file(resolve_path(params.db), type: "file", checkIfExists: true), "charon", "microbial", refs)
      query_ch = query_ch.mix(evaluate_charon.out.blast_queries)
      done_ch = done_ch.mix(reads_done(evaluate_charon.out))
    }

    if (params.charon_host_only){
      evaluate_charon_host_only(fastq_ch, file(resolve_path(params.db_host_only), type: "file", checkIfExists: true), "charon_host_only", "non_host", refs)
      query_ch = query_ch.mix(evaluate_charon_host_only.out.blast_queries)
      done_ch = done_ch.mix(reads_done(evaluate_charon_host_only.out))
    }

    if (params.deacon || params.deacon_sensitive){
      download_deacon_index()
    }

    if (params.deacon){
      evaluate_deacon(fastq_ch, download_deacon_index.out, "deacon", "", refs)
      query_ch = query_ch.mix(evaluate_deacon.out.blast_queries)
      done_ch = done_ch.mix(reads_done(evaluate_deacon.out))
    }

    if (params.deacon_sensitive){
      // deacon's high sensitivity settings
      evaluate_deacon_sensitive(fastq_ch, download_deacon_index.out, "deacon_sensitive", "-a 1 -r 0", refs)
      query_ch = query_ch.mix(evaluate_deacon_sensitive.out.blast_queries)
      done_ch = done_ch.mix(reads_done(evaluate_deacon_sensitive.out))
    }

    if (params.nohuman){
      evaluate_nohuman(fastq_ch, refs)
      query_ch = query_ch.mix(evaluate_nohuman.out.blast_queries)
      done_ch = done_ch.mix(reads_done(evaluate_nohuman.out))
    }

    if (params.evaluate_microbial) {
        blast(query_ch)
        blast_results = blast.out.results
    } else {
        blast_results = channel.empty()
    }

    eval_ch = channel.empty()
    if (params.charon)
        eval_ch = eval_ch.mix(build_eval_ch(evaluate_charon.out, "charon", blast_results))
    if (params.charon_host_only)
        eval_ch = eval_ch.mix(build_eval_ch(evaluate_charon_host_only.out, "charon_host_only", blast_results))
    if (params.deacon)
        eval_ch = eval_ch.mix(build_eval_ch(evaluate_deacon.out, "deacon", blast_results))
    if (params.deacon_sensitive)
        eval_ch = eval_ch.mix(build_eval_ch(evaluate_deacon_sensitive.out, "deacon_sensitive", blast_results))
    if (params.nohuman)
        eval_ch = eval_ch.mix(build_eval_ch(evaluate_nohuman.out, "nohuman", blast_results))

    // a process can only be invoked once per workflow, so all methods share this call
    contaminant_accs = params.contaminant_accs ? file(params.contaminant_accs, checkIfExists: true) : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    evaluate_summary(eval_ch, contaminant_accs)

    emit:
    done = done_ch
}

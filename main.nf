include { bam_to_fastq } from './modules/utils'
include { stage_fastq; build_mmi_index } from './modules/utils'
include { evaluate_dehosting } from './modules/evaluate_dehosting'
include { throttle_samples; release_finished_samples } from './modules/gate'

// Download remote inputs once; local files are used where they are.
workflow stage_remote {
    take:
    fastq_ch

    main:
    fastq_ch.branch{ id, f ->
                remote: f.toUri().scheme != 'file'
                local: true
            }
            .set{ branched }
    stage_fastq(branched.remote)
    out = branched.local.mix(stage_fastq.out)

    emit:
    out
}

// Slots left behind by a crashed or killed run would block this one, so clear them on start and on
// exit. This assumes slot_dir is used by one pipeline at a time (use a different --slot_dir per
// concurrent pipeline, or --reset_slots false to share a limit across pipelines).
workflow.onComplete {
    if (params.reset_slots) file(params.slot_dir).deleteDir()
}

workflow {
    if (params.reset_slots) file(params.slot_dir).deleteDir()

    // Check if required parameters are provided
    if (params.charon && !params.db) {
        exit 1, "Please provide --db (the charon index) when running Nextflow with --charon."
    }
    if (params.charon_host_only && !params.db_host_only) {
        exit 1, "Please provide --db_host_only (a charon index of the host only) when running Nextflow with --charon_host_only."
    }
    if (params.charon_host_only && !params.charon_host_only_container) {
        exit 1, "Please provide --charon_host_only_container, a charon image which supports a host-only index, when running with --charon_host_only."
    }

    if (params.read_type == "illumina") {
        exit 1, "I don't support paired end reads yet"
    }

    if (params.bam_csv){
        Channel.fromPath("${params.bam_csv}")
            .splitCsv( header: true )
            .map{ row -> [row.sample_id,file(row.filepath, type: "file", checkIfExists:true)]}
            .set{ bam_ch }
        throttle_samples(bam_ch)
        bam_to_fastq(throttle_samples.out)
        evaluate_dehosting(bam_to_fastq.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq_dir){
        fastq_ch = Channel.fromPath("${params.fastq_dir}/*.f*q*", type: "file", checkIfExists: true).map { [it.baseName.replace(".fq","").replace(".fastq",""), it] }
        throttle_samples(fastq_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq){
        unique_id = "${params.unique_id}"
        fastq = file(params.fastq, type: "file", checkIfExists:true)
        fastq_ch = Channel.from([[unique_id, fastq]])
        throttle_samples(fastq_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else {
        exit 1, "Please either provide a CSV with --bam_csv (specifying sample_id,filepath), a directory of fastq with --fastq_dir (each file for sample) or --fastq and --unique_id"
    }

    // free each sample's slot once every enabled method has finished with it
    n_methods = [params.charon, params.charon_host_only, params.deacon, params.deacon_sensitive, params.nohuman].count{ it }
    release_finished_samples(done_ch, n_methods)
}

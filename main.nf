include { bam_to_fastq } from './modules/utils'
include { stage_fastq; resolve_path; paired_input; csv_rows } from './modules/utils'
include { evaluate_dehosting } from './modules/evaluate_dehosting'
include { throttle_samples; release_finished_samples } from './modules/gate'

// Without a declared type, a command-line value is a string (`--charon false` would be true). Defaults are in
// conf/params.config.
params {
    paired: Boolean
    charon: Boolean
    charon_host_only: Boolean
    deacon: Boolean
    deacon_sensitive: Boolean
    nohuman: Boolean
    evaluate_host: Boolean
    evaluate_microbial: Boolean
    publish_reads: Boolean
    map_unclassified: Boolean
    blast_nonhuman: Boolean
    blast_cache: Boolean
    reset_slots: Boolean
    eager_cleanup: Boolean
    s3_anonymous: Boolean
    max_cpus: Integer
    host_eval_reads: Integer
    chunk_size: Integer
    max_concurrent_samples: Integer
    slot_timeout: Integer
    blast_max_targets: Integer
    blast_threads: Integer
    blast_chunk_reads: Integer?
    blast_max_reads: Integer?
    blast_min_cov: Float
    blast_min_pident: Float
    blast_bitscore_ratio: Float
}

workflow stage_remote {
    take:
    fastq_ch

    main:
    fastq_ch.branch{ id, f ->
                remote: (f instanceof List ? f : [f]).any{ it.toUri().scheme != 'file' }
                local: true
            }
            .set{ branched }
    stage_fastq(branched.remote)
    out = branched.local.mix(stage_fastq.out)

    emit:
    out
}

// Leftover slots from a killed run would block this one. Use a separate --slot_dir per concurrent pipeline.
workflow {
    main:
    if (params.reset_slots) file(params.slot_dir).deleteDir()

    if (params.charon && !params.db) {
        error "Please provide --db (the charon index) when running Nextflow with --charon."
    }
    if (params.charon_host_only && !params.db_host_only) {
        error "Please provide --db_host_only (a charon index of the host only) when running Nextflow with --charon_host_only."
    }

    // FASTQ pairs are recognised by having a second file; a BAM/CRAM of pairs (--bam_csv) needs --paired
    if (params.paired && !params.bam_csv) {
        error "--paired is only for --bam_csv: FASTQ pairs are recognised by their second file (--fastq2, or fastq_2 in --fastq_csv)."
    }
    if (paired_input() && params.read_type != "illumina") {
        log.info "Paired reads are mapped as short reads (the sr preset), whatever --read_type is."
    }

    if (params.bam_csv){
        channel.fromPath("${params.bam_csv}")
            .splitCsv( header: true, quote: '"' )
            .map{ row -> [row.sample_id,file(row.filepath, type: "file", checkIfExists:true)]}
            .set{ bam_ch }
        cram_ref = params.cram_ref ? file(resolve_path(params.cram_ref), type: "file", checkIfExists: true) : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
        throttle_samples(bam_ch)
        bam_to_fastq(throttle_samples.out, cram_ref)
        evaluate_dehosting(bam_to_fastq.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq_csv){
        channel.fromList(csv_rows())
            .map{ row ->
                def r1 = file(row.fastq_1, type: "file", checkIfExists: true)
                [row.sample_id, row.fastq_2 ? [r1, file(row.fastq_2, type: "file", checkIfExists: true)] : r1]
            }
            .set{ fastq_csv_ch }
        throttle_samples(fastq_csv_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq_dir){
        fastq_ch = channel.fromPath("${params.fastq_dir}/*.f*q*", type: "file", checkIfExists: true).map { [it.baseName.replace(".fq","").replace(".fastq",""), it] }
        throttle_samples(fastq_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq){
        fastq = file(params.fastq, type: "file", checkIfExists:true)
        // a String, as a GString key would not match in joins
        unique_id = params.unique_id ? params.unique_id.toString() : fastq.name.replaceAll(/\.(fq|fastq)(\.gz)?$/, '')
        reads = params.fastq2 ? [fastq, file(params.fastq2, type: "file", checkIfExists:true)] : fastq
        fastq_ch = channel.from([[unique_id, reads]])
        throttle_samples(fastq_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else {
        error "Please either provide a CSV with --bam_csv (specifying sample_id,filepath), a CSV of FASTQ with --fastq_csv (sample_id,fastq_1[,fastq_2]), a directory of fastq with --fastq_dir (each file for sample) or --fastq (and --fastq2), with --unique_id to name the sample"
    }

    n_methods = [params.charon, params.charon_host_only, params.deacon, params.deacon_sensitive, params.nohuman].count{ it }
    release_finished_samples(done_ch, n_methods)

    onComplete:
    if (params.reset_slots) file(params.slot_dir).deleteDir()
}

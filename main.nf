include { bam_to_fastq } from './modules/utils'
include { stage_fastq; resolve_path; paired_input } from './modules/utils'
include { evaluate_dehosting } from './modules/evaluate_dehosting'
include { throttle_samples; release_finished_samples } from './modules/gate'

// Types of the boolean and numeric parameters. Nextflow 26 passes command-line values to the script as strings
// unless the parameter is declared here, so `--charon false` would count as true and `--host_eval_reads 10` would
// not be a number. Defaults are in conf/params.config.
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
    charon_chunk_size: Integer
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

// Download remote inputs once; local files are used where they are.
workflow stage_remote {
    take:
    fastq_ch

    main:
    fastq_ch.branch{ id, f ->
                // f is one file, or a list of two for a pair
                remote: (f instanceof List ? f : [f]).any{ it.toUri().scheme != 'file' }
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
workflow {
    main:
    if (params.reset_slots) file(params.slot_dir).deleteDir()

    // Check if required parameters are provided
    if (params.charon && !params.db) {
        exit 1, "Please provide --db (the charon index) when running Nextflow with --charon."
    }
    if (params.charon_host_only && !params.db_host_only) {
        exit 1, "Please provide --db_host_only (a charon index of the host only) when running Nextflow with --charon_host_only."
    }

    // Paired reads: FASTQ pairs (--fastq_pairs_csv, or --fastq and --fastq2) are paired by themselves; a BAM/CRAM
    // of pairs (--bam_csv) needs --paired, because the file does not say whether it holds pairs.
    if (params.paired && (params.fastq_dir || (params.fastq && !params.fastq2))) {
        exit 1, "--paired needs --bam_csv, --fastq_pairs_csv, or --fastq with --fastq2 (--fastq_dir is not supported for pairs)."
    }
    if (paired_input() && params.read_type != "illumina") {
        log.info "Paired reads are mapped as short reads (the sr preset), whatever --read_type is."
    }

    if (params.bam_csv){
        Channel.fromPath("${params.bam_csv}")
            .splitCsv( header: true, quote: '"' )
            .map{ row -> [row.sample_id,file(row.filepath, type: "file", checkIfExists:true)]}
            .set{ bam_ch }
        // a CRAM needs the reference it was aligned to; a BAM (or an unaligned file) does not
        cram_ref = params.cram_ref ? file(resolve_path(params.cram_ref), type: "file", checkIfExists: true) : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
        throttle_samples(bam_ch)
        bam_to_fastq(throttle_samples.out, cram_ref)
        evaluate_dehosting(bam_to_fastq.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq_pairs_csv){
        Channel.fromPath("${params.fastq_pairs_csv}")
            .splitCsv( header: true, quote: '"' )
            .map{ row -> [row.sample_id, [file(row.fastq_1, type: "file", checkIfExists:true), file(row.fastq_2, type: "file", checkIfExists:true)]]}
            .set{ pairs_ch }
        throttle_samples(pairs_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq_dir){
        fastq_ch = Channel.fromPath("${params.fastq_dir}/*.f*q*", type: "file", checkIfExists: true).map { [it.baseName.replace(".fq","").replace(".fastq",""), it] }
        throttle_samples(fastq_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else if (params.fastq){
        fastq = file(params.fastq, type: "file", checkIfExists:true)
        // a plain String: a GString key would not match the String keys used in joins further down.
        // Without --unique_id, the name of the (first) FASTQ file without .fq/.fastq(.gz)
        unique_id = params.unique_id ? params.unique_id.toString() : fastq.name.replaceAll(/\.(fq|fastq)(\.gz)?$/, '')
        reads = params.fastq2 ? [fastq, file(params.fastq2, type: "file", checkIfExists:true)] : fastq
        fastq_ch = Channel.from([[unique_id, reads]])
        throttle_samples(fastq_ch)
        stage_remote(throttle_samples.out)
        evaluate_dehosting(stage_remote.out)
        done_ch = evaluate_dehosting.out.done
    } else {
        exit 1, "Please either provide a CSV with --bam_csv (specifying sample_id,filepath), a CSV of FASTQ pairs with --fastq_pairs_csv (sample_id,fastq_1,fastq_2), a directory of fastq with --fastq_dir (each file for sample) or --fastq (and --fastq2), with --unique_id to name the sample"
    }

    // free each sample's slot once every enabled method has classified and mapped it (BLAST and the summaries come
    // later, for all samples together)
    n_methods = [params.charon, params.charon_host_only, params.deacon, params.deacon_sensitive, params.nohuman].count{ it }
    release_finished_samples(done_ch, n_methods)

    onComplete:
    if (params.reset_slots) file(params.slot_dir).deleteDir()
}

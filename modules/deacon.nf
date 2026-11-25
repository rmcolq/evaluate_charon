#!/usr/bin/env nextflow
include { minimap2_microbial; minimap2_host; extract_blast_queries; publish_reads; is_paired; reads_prefix; reads_glob; first_file } from '../modules/utils'

process download_deacon_index {
    label "process_single"
    storeDir "${params.store_dir}/deacon/"
    container 'community.wave.seqera.io/library/deacon:0.20.0--903cc7096f1e9938'
    maxForks 1

    output:
        path("*.idx")

    script:
    """
    export SSL_CERT_FILE=/opt/conda/ssl/cacert.pem
    deacon index fetch ${params.deacon_index} -k 31 -w 15
    """
}
process run_deacon {

    label "process_medium"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/deacon:0.20.0--903cc7096f1e9938'
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
    """
    deacon filter ${extra_args} ${deacon_index} ${fastq} ${deacon_output_args(name, unique_id, 'human', fastq)}
    deacon filter -d ${extra_args} ${deacon_index} ${fastq} ${deacon_output_args(name, unique_id, 'microbial', fastq)}
    """
}

// deacon options naming the output file(s) of a category: -o, plus -O for the second mate of a pair
def deacon_output_args(name, unique_id, category, reads) {
    return is_paired(reads) ?
        "-o ${name}_${unique_id}_${category}_1.fq.gz -O ${name}_${unique_id}_${category}_2.fq.gz" :
        "-o ${name}_${unique_id}_${category}.fq.gz"
}

// Output file(s) of a category: one file, or _1 and _2 for a pair
def deacon_out(name, unique_id, category, reads) {
    return is_paired(reads) ? "${name}_${unique_id}_${category}_?.fq.gz" : "${name}_${unique_id}_${category}.fq.gz"
}

process collect_classifications {

    label "process_low"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/deacon:0.20.0--903cc7096f1e9938'

    publishDir path: { "${params.outdir}/${unique_id}/" }, mode: 'copy', pattern: "*.out"

    input:
    tuple val(unique_id), val(name), path(microbial_fastq), path(host_fastq)

    output:
    tuple val(unique_id), val(name), path("${unique_id}.${name}.out"),  emit: result

    script:
    """
    # one line per read, or per pair (the id taken from R1, without /1 or /2)
    echo -e "read_id\tclassification" > "${unique_id}.${name}.out"
    zcat -f ${first_file(microbial_fastq)} | awk 'NR%4==1 {id=substr(\$1,2); sub(/\\/[12]\$/, "", id); print id "\tmicrobial"}' >> "${unique_id}.${name}.out"
    zcat -f ${first_file(host_fastq)} | awk 'NR%4==1 {id=substr(\$1,2); sub(/\\/[12]\$/, "", id); print id "\thuman"}' >> "${unique_id}.${name}.out"
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
        refs
    main:


    run_deacon(fastq_ch, deacon_index, name, extra_args)
    if ( params.publish_reads ) {
        publish_reads(run_deacon.out.microbial_fastq.map{ id, method, reads -> [id, method, "microbial", reads] })
    }
    collect_classifications(run_deacon.out.combined)
    minimap2_microbial(run_deacon.out.microbial_fastq, refs)
    minimap2_host(run_deacon.out.human_fastq, refs)

    // BLAST is run once for all samples and methods (modules/blast.nf, called from modules/evaluate_dehosting.nf), so only the
    // queries are extracted here
    if ( params.evaluate_microbial ){
        extract_blast_queries(minimap2_microbial.out.chunk_sam_ch)
        query_ch = extract_blast_queries.out.queries
    } else {
        query_ch = channel.empty()
    }

    emit:
        report = collect_classifications.out.result
        microbial_sam = minimap2_microbial.out.sam_ch
        host_sam = minimap2_host.out
        blast_queries = query_ch
        // samples with no microbial reads have no queries
        no_reads = minimap2_microbial.out.no_reads
}

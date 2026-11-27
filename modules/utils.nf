#!/usr/bin/env nextflow

// A path parameter: relative to the launch directory, else to the pipeline (resources/). A URL is downloaded once into
// <store_dir>/downloads.
def resolve_path(String path) {
    if (path.contains('://')) {
        def name = path.replaceAll(/\?.*$/, '').replaceAll(/\/content$/, '').tokenize('/').last()
        def dest = file("${params.store_dir}/downloads/${path.md5().take(8)}/${name}")
        if (!dest.exists()) {
            log.info "Downloading ${path} to ${dest}"
            def tmp = file("${dest}.part")
            tmp.parent.mkdirs()
            file(path).copyTo(tmp)
            tmp.moveTo(dest)
        }
        return dest
    }
    def local = file(path)
    return (path.startsWith('/') || local.exists()) ? local : file("${projectDir}/${path}")
}

// A path option which the tool needs uncompressed (the Charon indexes are deposited gzipped): a .gz file is
// decompressed once into <store_dir>/unzipped and that copy is used
def resolve_unzipped(String path) {
    def f = resolve_path(path)
    if (!f.name.endsWith(".gz"))
        return f
    def dest = file("${params.store_dir}/unzipped/${f.toString().md5().take(8)}/${f.name.replaceFirst(/\.gz$/, '')}")
    if (!dest.exists()) {
        log.info "Decompressing ${f} to ${dest}"
        def tmp = file("${dest}.part")
        tmp.parent.mkdirs()
        new java.util.zip.GZIPInputStream(f.newInputStream(), 1 << 20).withCloseable { input ->
            java.nio.file.Files.copy(input, tmp, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
        }
        tmp.moveTo(dest)
    }
    return dest
}

// All rows must be single-end or all paired: a run maps with one minimap2 preset
def csv_rows() {
    def rows = file(params.fastq_csv, checkIfExists: true).splitCsv(header: true, quote: '"')
    def paired = rows.collect{ row -> row.fastq_2 ? true : false }.unique()
    if (paired.size() > 1)
        error "--fastq_csv mixes single-end and paired samples: run them separately"
    return rows
}

def paired_input() {
    return params.paired || params.fastq2 || (params.fastq_csv && csv_rows().any{ row -> row.fastq_2 })
}

def minimap2_preset() {
    return (paired_input() || params.read_type == "illumina") ? "sr" : "map-ont"
}

// Nanopore reads have a far higher error rate than Illumina: leaked human reads align at ~75% identity
// (Hall et al. artificial ONT), so the BLAST search and leak tier 3 use looser defaults for them.
def noisy_reads() {
    return minimap2_preset() == "map-ont"
}

def blast_min_pident() {
    return params.blast_min_pident != null ? params.blast_min_pident : (noisy_reads() ? 70 : 90)
}

def blast_min_cov() {
    return params.blast_min_cov != null ? params.blast_min_cov : (noisy_reads() ? 50 : 80)
}

// dc-megablast seeds on discontiguous words, so it finds hits megablast misses in error-rich reads
def blast_task() {
    return params.blast_task ?: (noisy_reads() ? "dc-megablast" : "megablast")
}

// identity filter applied during the search: no stricter than tier 3 needs
def blast_search_pident() {
    return Math.min(85, blast_min_pident() as double)
}

// Reads travel as one file, or a list of two (R1, R2) for a pair
def is_paired(reads) {
    return reads instanceof List
}

def reads_prefix(unique_id, reads) {
    return is_paired(reads) ? "${unique_id}" : reads.baseName
}

def first_file(reads) {
    return is_paired(reads) ? reads[0] : reads
}

def reads_glob(prefix, category, reads) {
    return is_paired(reads) ? "${prefix}_${category}_?.f*q.gz" : "${prefix}_${category}.f*q.gz"
}

def reads_name(reads) {
    return is_paired(reads) ? reads[0].baseName : reads.baseName
}

process bam_to_fastq {
    label "process_medium"
    tag "${unique_id}"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"
    maxForks 2

    input:
    tuple val(unique_id), path(bam)
    path(cram_ref)

    output:
    tuple val(unique_id), path(params.paired ? "${unique_id}_R*.fastq.gz" : "${bam.baseName}.fastq.gz")

    script:
    def ref = cram_ref.name == "NO_FILE" ? "" : "--reference ${cram_ref}"
    if (params.paired) {
        """
        # primary alignments only, grouped by read name so that mates are adjacent
        samtools view -u ${ref} -F 0x900 ${bam} \
          | samtools collate -O -u - collate_tmp \
          | samtools fastq -n -1 "${unique_id}_R1.fastq.gz" -2 "${unique_id}_R2.fastq.gz" -0 /dev/null -s /dev/null -
        """
    } else {
        """
        # host evaluation and charon's model use the first reads, which in a sorted file all come from chr1
        sorted=\$(samtools view -H ${ref} ${bam} | grep -c "SO:coordinate" || true)
        if [ "\$sorted" -gt 0 ]; then
          samtools collate -O -u -@ ${task.cpus} ${ref} ${bam} collate_tmp | samtools fastq -
        else
          samtools fastq ${ref} ${bam}
        fi | bgzip -@ ${task.cpus} > "${bam.baseName}.fastq.gz"
        """
    }
}

process stage_fastq {
    // downloads a remote input once, rather than once per classifier
    label "process_single"
    tag "${unique_id}"
    maxForks 2

    input:
    tuple val(unique_id), path(fastq)

    output:
    tuple val(unique_id), path("staged_*")

    script:
    """
    for f in ${fastq}; do mv "\$f" "staged_\$f"; done
    """
}

process shuffle_fastq {
    // Charon fits its model on the first reads, so a file grouped by source genome (as the Hall et al. mixtures
    // are) is shuffled first. In memory (seqkit's two-pass mode is FASTA only); the same seed shuffles both files
    // of a pair the same way, so mates stay in step. seqkit is pinned, as another version may give another order.
    label "process_medium_plus_mem"
    label "error_retry"
    container "quay.io/biocontainers/seqkit:2.14.0--hb192632_0"
    tag "${unique_id}"

    input:
    tuple val(unique_id), path(fastq)

    output:
    tuple val(unique_id), path("shuffled_*")

    script:
    def files = fastq instanceof List ? fastq : [fastq]
    def cmds = files.collect{ f ->
        "seqkit shuffle -s ${params.shuffle_seed} -j ${task.cpus} \"${f}\" -o \"shuffled_${f.name.replaceFirst(/^staged_/, '')}\""
    }.join("\n")
    """
    ${cmds}
    """
}

// Built once into the store dir. The reference size is in the name, so a changed reference is indexed again.
process minimap2_index {
    label "process_medium"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"
    storeDir "${params.store_dir}/minimap2/"
    maxForks 1

    input:
        path refs
        val index_name
    output:
        path(index_name)

    script:
        """
        minimap2 -x ${minimap2_preset()} -t ${task.cpus} -d "${index_name}" ${refs}
        """
}

workflow mapping_index {
    main:
    if ( params.evaluate_host || params.evaluate_microbial ) {
        refs = file(resolve_path(params.refs), type: "file", checkIfExists: true)
        index_name = "${refs.simpleName}.${refs.size()}.${minimap2_preset()}.mmi"
        log.info "minimap2 index: ${file(params.store_dir).toAbsolutePath()}/minimap2/${index_name} (built if missing)"
        minimap2_index(refs, index_name)
        index_ch = minimap2_index.out
    } else {
        index_ch = channel.value(file("${projectDir}/assets/NO_FILE", checkIfExists: true))
    }

    emit:
    index = index_ch
}

process minimap2_microbial_chunk {

    label "process_medium"
    tag "${unique_id}"

    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), val(method), path("microbial.${reads_name(fastq)}.mmp.sam")
    script:
        if ( params.evaluate_microbial ) {
            preset = minimap2_preset()
            """
            minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > microbial.${reads_name(fastq)}.mmp.sam
            """
        } else {
            """
            touch "microbial.${reads_name(fastq)}.mmp.sam"
            """
        }

}

// R1 and R2 are split the same way, so chunk i of each holds the same pairs
def split_reads(reads) {
    if (!is_paired(reads))
        return reads.splitFastq(by: params.chunk_size, file: true)
    def r1 = reads[0].splitFastq(by: params.chunk_size, file: true)
    def r2 = reads[1].splitFastq(by: params.chunk_size, file: true)
    return [r1, r2].transpose()
}

workflow minimap2_microbial {
    take:
    fastq_ch
    refs

    main:
    def empty_sam = file("${projectDir}/resources/empty.sam", type: "file", checkIfExists: true)

    // a sample with no microbial reads has no chunks, so it gets an empty result or it would never finish
    fastq_ch.map{unique_id, method, fastq -> [unique_id, method, split_reads(fastq)]}
            .branch{ unique_id, method, chunks ->
                with_reads: chunks.size() > 0
                no_reads: true
            }
            .set{ split_ch }
    // groupKey lets groupTuple emit a sample when its own chunks are done; waiting for the channel to close would
    // deadlock with the sample gate
    split_ch.with_reads
        .map{ unique_id, method, chunks -> [groupKey(unique_id, chunks.size()), method, chunks] }
        .transpose()
        .set{chunked_fastq_ch}
    minimap2_microbial_chunk(chunked_fastq_ch, refs)
    minimap2_microbial_chunk.out.set{ chunk_sam_ch }
    minimap2_microbial_chunk.out
        .groupTuple()
        .map{ unique_id, methods, sams -> [unique_id, methods[0], sams] }
        .set{collected_sam}
    cat_all_microbial_sam_files(collected_sam)
    cat_all_microbial_sam_files.out
        .map{ unique_id, sam -> [unique_id.toString(), sam] }
        .mix( split_ch.no_reads.map{ unique_id, method, chunks -> [unique_id, empty_sam] } )
        .set{ sam_ch }
    split_ch.no_reads.map{ unique_id, method, chunks -> [unique_id, method] }.set{ no_reads }

    emit:
    chunk_sam_ch
    sam_ch
    no_reads
}

process minimap2_host {

    label "process_medium"
    tag "${unique_id}"

    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    publishDir path: { "${params.outdir}/${unique_id}/mapped/" }, mode: 'copy', pattern: "*.sam"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), path("${unique_id}.${method}.host.sam"), path("${unique_id}.${method}.host_sampled.txt")
    script:
        if ( params.evaluate_host == true ) {
            preset = minimap2_preset()
            inputs = is_paired(fastq) ? "small_1.fq small_2.fq" : "small.fq"
            """
            # head closes the pipe once it has enough reads, which gives zcat SIGPIPE (141) under pipefail
            i=0
            for f in ${fastq}; do
              i=\$((i+1))
              zcat -f \$f | head -n ${params.host_eval_reads * 4} > small_\$i.fq || [ \$? -eq 141 ]
            done
            ${is_paired(fastq) ? "" : "mv small_1.fq small.fq"}
            # the sampled reads, so that reads not sampled are not counted as unmapped
            awk 'NR%4==1 {print substr(\$1,2)}' ${is_paired(fastq) ? "small_1.fq" : "small.fq"} | sed 's#/[12]\$##' > "${unique_id}.${method}.host_sampled.txt"
            minimap2 -ax ${preset} ${refs} ${inputs} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > "${unique_id}.${method}.host.sam"
            """
        } else {
            """
            touch "${unique_id}.${method}.host.sam" "${unique_id}.${method}.host_sampled.txt"
            """
        }
}

process extract_microbial_host_hits {

    label "process_single"
    tag "${unique_id}"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"

    input:
    tuple val(unique_id), val(method), path(sam_file)
    path ref_bed

    output:
    // unique names, as the files of all samples and methods are staged together for BLAST
    tuple val(unique_id), val(method), path("${unique_id}.${method}.${sam_file.baseName}.query.fasta")

    script:
    """
    # both mates of a pair are kept, as /1 and /2; evaluate.py takes the more human one
    samtools view -u -L ${ref_bed} ${sam_file} | samtools fasta - > "${unique_id}.${method}.${sam_file.baseName}.query.fasta"
    """
}

// Maps all the reads charon leaves unclassified, not a sample
process minimap2_unclassified {

    label "process_medium"
    label "error_ignore"
    tag "${unique_id}"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"
    publishDir path: { "${params.outdir}/${unique_id}/mapped/" }, mode: 'copy', pattern: "*.sam"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), val(method), path("${unique_id}.${method}.unclassified.sam")
    script:
        """
        minimap2 -ax ${minimap2_preset()} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > "${unique_id}.${method}.unclassified.sam"
        """
}

process alignment_table {

    label "process_low"
    label "error_ignore"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    publishDir path: { "${params.outdir}/${unique_id}/" }, mode: 'copy'

    input:
        tuple val(unique_id), val(method), path(sam)
    output:
        path "${unique_id}_${method}_unclassified_alignments.csv.gz"
    script:
        """
        alignment_table.py --sam ${sam} --out "${unique_id}_${method}_unclassified_alignments.csv.gz"
        """
}

process publish_reads {

    label "process_single"
    label "error_ignore"
    tag "${unique_id}"
    publishDir path: { "${params.outdir}/${unique_id}/reads/" }, mode: 'copy'

    input:
    tuple val(unique_id), val(method), val(category), path(reads, stageAs: "in/*")

    output:
    path "${unique_id}_${method}_${category}*.f*q.gz"

    script:
    def files = reads instanceof List ? reads : [reads]
    def copies = files.withIndex().collect{ f, i ->
        def suffix = files.size() > 1 ? "_${i + 1}" : ""
        def out = "${unique_id}_${method}_${category}${suffix}.fastq.gz"
        "cp -L \"${f}\" tmp.gz && (gzip -t tmp.gz 2>/dev/null && mv tmp.gz ${out} || gzip -c \"${f}\" > ${out})"
    }
    """
    ${copies.join("\n    ")}
    """
}

process cat_all_microbial_sam_files {

    label "process_low"
    tag "${unique_id}"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"
    publishDir path: { "${params.outdir}/${unique_id}/mapped/" }, mode: 'copy', pattern: "*.sam"

    input:
    tuple val(unique_id), val(method), path(sam_files)

    output:
    tuple val(unique_id), path("${unique_id}.${method}.microbial.sam")

    script:
    """
    cat \$(ls *.mmp.sam | head -n1) > "${unique_id}.${method}.microbial.sam"
    for sam in \$(ls *.mmp.sam | tail -n+2)
      do
        grep -v '^@' \$sam >> "${unique_id}.${method}.microbial.sam" || true
      done
    """
}


workflow extract_blast_queries {
    take:
        chunk_sam_ch

    main:
        ref_bed = file(resolve_path(params.ref_bed), type: "file", checkIfExists:true)
        extract_microbial_host_hits(chunk_sam_ch, ref_bed)
        // the sample id arrives as a groupKey from the chunking; later joins need plain String keys
        extract_microbial_host_hits.out
            .map{ unique_id, method, query -> [unique_id.toString(), method, query] }
            .set{ queries }

    emit:
        queries
}

process evaluate_summary {

    label "process_low"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    publishDir path: { "${params.outdir}/${unique_id}/" }, mode: 'copy'

    input:
    tuple val(unique_id), val(classifier), path(report), path(host_sam), path(host_sampled_ids), path(microbial_sam), path(blast_result)
    path(contaminant_accs)

    output:
    tuple val(unique_id), val(classifier), emit: done
    path "${unique_id}*_summary.csv", emit: summary
    path "${unique_id}*_full.csv", emit: full
    path "${unique_id}*_calls.*", emit: calls
    path "${unique_id}*_leaks.*", emit: leaks
    path "${unique_id}*_data.csv", emit: data
    path "${unique_id}*_taxa.csv", emit: taxa, optional:true
    path "${unique_id}*_accs.csv", emit: accs, optional:true

    script:
    """
    evaluate.py \
      -i ${report} \
      --microbial_sam ${microbial_sam} \
      --host_sam ${host_sam} \
      --host_sampled_ids ${host_sampled_ids} \
      ${params.evaluate_host ? "" : "--no_host_mapping"} \
      ${params.evaluate_microbial ? "" : "--no_microbial_mapping"} \
      --blast_result ${blast_result} \
      --classifier ${classifier} \
      --sample_id ${unique_id} \
      --contaminant_accs ${contaminant_accs} \
      --min_blast_cov ${blast_min_cov()} \
      --min_blast_pident ${blast_min_pident()} \
      --min_aligned_prop ${params.leak_min_aligned_prop} \
      --min_mapq ${params.leak_min_mapq} \
      --bitscore_ratio ${params.blast_bitscore_ratio} \
      -p "${unique_id}_${classifier}"
    """
}
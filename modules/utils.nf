#!/usr/bin/env nextflow

// A path given as a parameter: absolute, or relative to the pipeline directory.
def resolve_path(String path) {
    // test the text, not file(path), which would already resolve a relative path against the launch directory
    return (path.startsWith('/') || path.contains('://')) ? file(path) : file("${projectDir}/${path}")
}

// Paired input: FASTQ pairs are paired by their shape; a BAM/CRAM only when --paired says so.
def paired_input() {
    return params.paired || params.fastq2 || params.fastq_pairs_csv
}

def minimap2_preset() {
    return (paired_input() || params.read_type == "illumina") ? "sr" : "map-ont"
}

// Reads travel through the pipeline as one channel item: a single file, or a list of two files (R1, R2)
// for paired reads. Each process works out which it has from the shape of its input.
def is_paired(reads) {
    return reads instanceof List
}

// Prefix of charon's output files: the sample for a pair, the file name for single-end reads
def reads_prefix(unique_id, reads) {
    return is_paired(reads) ? "${unique_id}" : reads.baseName
}

// The first (or only) file of a set of reads
def first_file(reads) {
    return is_paired(reads) ? reads[0] : reads
}

// Glob of the extracted reads files for a category: one file, or _1 and _2 for a pair
def reads_glob(prefix, category, reads) {
    return is_paired(reads) ? "${prefix}_${category}_?.f*q.gz" : "${prefix}_${category}.f*q.gz"
}

// Name of the first (or only) reads file without its extensions
def reads_name(reads) {
    return is_paired(reads) ? reads[0].baseName : reads.baseName
}

process bam_to_fastq {
    // Reads from a BAM or CRAM. With --paired the file holds read pairs, which are written as R1/R2 FASTQ
    // (reads without a mate are dropped). A CRAM needs its reference, --cram_ref.
    label "process_medium"
    tag "${unique_id}"
    conda "bioconda::samtools=1.21"
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
        # A coordinate-sorted file is first grouped by read name, which puts the reads in an order unrelated to their
        # position: host evaluation maps the first reads, which would otherwise all come from the start of chr1.
        # samtools fastq leaves out secondary and supplementary alignments.
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
    // Remote (e.g. s3://) inputs are downloaded by Nextflow into this task's work dir once,
    // rather than separately by every classifier. Eager cleanup removes it after last use.
    // Takes one file or a pair; a pair keeps its order (R1 then R2) because the names sort that way.
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

// The minimap2 index of the reference is built once and kept in the store dir, so every mapping task loads it
// instead of rebuilding it from the FASTA. Natively a human-sized index takes about 2 minutes and 3.5 GB, but in an
// emulated or small-memory container (e.g. amd64 on Apple Silicon) it can be very slow or be killed, so on a laptop
// build it natively and copy it to the path logged at the start of the run.
// The size of the reference is in the name, so a changed reference is indexed again rather than using a stale index.
process minimap2_index {
    label "process_medium"
    conda "bioconda::minimap2=2.28"
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

// Indexes the reference once for all methods. Nothing is built when no mapping is evaluated.
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

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), val(method), path("microbial.${reads_name(fastq)}.mmp.sam")
    script:
        // fastq is one file, or R1 and R2 of a pair (mapped together, so each pair is mapped once)
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

// Split a reads file, or a pair, into chunks of params.chunk_size reads. For a pair, R1 and R2 are split
// the same way, so chunk i of each holds the same read pairs.
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

    // A sample with no microbial reads (e.g. a human sample) has no chunks, so nothing would reach the
    // summary and the sample would never finish. Give it an empty result instead.
    // Chunks are lists of one file, or of two files (R1, R2) with the same reads, for paired reads
    fastq_ch.map{unique_id, method, fastq -> [unique_id, method, split_reads(fastq)]}
            .branch{ unique_id, method, chunks ->
                with_reads: chunks.size() > 0
                no_reads: true
            }
            .set{ split_ch }
    // Each chunk carries its sample's chunk count (groupKey) so that groupTuple can emit a sample as soon as
    // its own chunks are done. Without it, groupTuple waits for the whole channel to close, which never
    // happens while later samples are held back by the sample gate, so the sample never finishes.
    // This workflow is called once per method, so grouping on the sample alone is enough.
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

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    publishDir path: { "${params.outdir}/${unique_id}/mapped/" }, mode: 'copy', pattern: "*.sam"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), path("${unique_id}.${method}.host.sam"), path("${unique_id}.${method}.host_sampled.txt")
    script:
        // fastq is one file, or R1 and R2 of a pair. host_eval_reads counts reads, or pairs for paired reads.
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
            # the reads which are mapped, so that the others are not counted as unmapped
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
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"

    input:
    tuple val(unique_id), val(method), path(sam_file)
    path ref_bed

    output:
    // named for the sample, method and chunk, so the files of all of them can be staged together for pooling
    tuple val(unique_id), val(method), path("${unique_id}.${method}.${sam_file.baseName}.query.fasta")

    script:
    """
    # -L streams through the alignments, keeping those overlapping the human reference, so no sorting or index is needed.
    # Both mates of a pair are kept, named /1 and /2 (single-end reads keep their name); evaluate.py takes the more
    # human mate for the pair.
    samtools view -u -L ${ref_bed} ${sam_file} | samtools fasta - > "${unique_id}.${method}.${sam_file.baseName}.query.fasta"
    """
}

// Map a set of reads which were not called by a method (e.g. the reads charon leaves unclassified) to the reference
// in full, not a sample, and publish the alignments. Reads which do not map are not in the SAM.
process minimap2_unclassified {

    label "process_medium"
    label "error_ignore"
    tag "${unique_id}"
    conda "bioconda::minimap2=2.28"
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

// One row per mapped read, with the same alignment columns as the full tables
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

// Keep a set of reads (one file, or R1 and R2 of a pair) in <outdir>/<sample>/reads, named for the sample, method and
// category. The reads are copied, because the originals are still needed by the next steps.
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


// The reads of one method that were called microbial but map to the human reference, ready for BLAST
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
      --min_blast_cov ${params.blast_min_cov} \
      --min_blast_pident ${params.blast_min_pident} \
      --bitscore_ratio ${params.blast_bitscore_ratio} \
      -p "${unique_id}_${classifier}"
    """
}
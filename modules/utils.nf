#!/usr/bin/env nextflow

// A path given as a parameter: absolute, or relative to the pipeline directory.
def resolve_path(String path) {
    // test the text, not file(path), which would already resolve a relative path against the launch directory
    return (path.startsWith('/') || path.contains('://')) ? file(path) : file("${projectDir}/${path}")
}

// The reference minimap2 maps to: a prebuilt index (--refs_mmi) if given, otherwise the FASTA (--refs).
// A prebuilt index avoids rebuilding the index in every mapping task, which needs several times more
// memory than mapping (about 12 GB for a human reference).
def mapping_reference() {
    if (params.refs_mmi)
        return file(resolve_path(params.refs_mmi), type: "file", checkIfExists: true)
    return file(resolve_path(params.refs), type: "file", checkIfExists: true)
}

def minimap2_preset() {
    return params.read_type == "illumina" ? "sr" : "map-ont"
}

process bam_to_fastq {
    label "process_medium"
    tag "${unique_id}"
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"
    maxForks 2

    input:
    tuple val(unique_id), path(bam)

    output:
    tuple val(unique_id), path("${bam.baseName}.fastq")

    script:
    """
    samtools fastq ${bam} > "${bam.baseName}.fastq"
    """
}

process stage_fastq {
    // Remote (e.g. s3://) inputs are downloaded by Nextflow into this task's work dir once,
    // rather than separately by every classifier. Eager cleanup removes it after last use.
    label "process_single"
    tag "${unique_id}"
    maxForks 2

    input:
    tuple val(unique_id), path(fastq)

    output:
    tuple val(unique_id), path("staged_${fastq}")

    script:
    """
    mv ${fastq} staged_${fastq}
    """
}

process build_mmi {
    // Builds the minimap2 index once (needs far more memory than mapping), so it can be reused with --refs_mmi.
    // Run with: nextflow run main.nf -entry build_mmi_index --refs <fasta> -profile ...
    label "process_high_memory"
    label "process_high"
    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    publishDir "${params.outdir}/", mode: 'copy', pattern: "*.mmi"

    input:
        path refs
    output:
        path("${refs.simpleName}.${minimap2_preset()}.mmi")

    script:
        preset = minimap2_preset()
        """
        minimap2 -x ${preset} -t ${task.cpus} -d "${refs.simpleName}.${preset}.mmi" ${refs}
        """
}

workflow build_mmi_index {
    main:
    build_mmi(file(resolve_path(params.refs), type: "file", checkIfExists: true))
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
        tuple val(unique_id), val(method), path("microbial.${fastq.baseName}.mmp.sam")
    script:
        if ( params.evaluate_microbial ) {
            preset = minimap2_preset()
            """
            minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > microbial.${fastq.baseName}.mmp.sam
            """
        } else {
            """
            touch "microbial.${fastq.baseName}.mmp.sam"
            """
        }
        
}

workflow minimap2_microbial {
    take:
    fastq_ch
    refs

    main:
    def empty_sam = file("${projectDir}/resources/empty.sam", type: "file", checkIfExists: true)

    // A sample with no microbial reads (e.g. a human sample) has no chunks, so nothing would reach the
    // summary and the sample would never finish. Give it an empty result instead.
    fastq_ch.map{unique_id, method, fastq -> [unique_id, method, fastq.splitFastq(by: params.chunk_size, file:true)]}
            .branch{ unique_id, method, chunks ->
                with_reads: chunks.size() > 0
                no_reads: true
            }
            .set{ split_ch }
    split_ch.with_reads.transpose().set{chunked_fastq_ch}
    minimap2_microbial_chunk(chunked_fastq_ch, refs)
    minimap2_microbial_chunk.out.set{ chunk_sam_ch }
    minimap2_microbial_chunk.out.groupTuple(by: [0,1]).set{collected_sam}
    cat_all_microbial_sam_files(collected_sam)
    cat_all_microbial_sam_files.out
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

    publishDir "${params.outdir}/intermediate/", mode: 'copy', pattern: "*.sam"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), path("${unique_id}.${method}.host.sam"), path("${unique_id}.${method}.host_sampled.txt")
    script:
        if ( params.evaluate_host == true ) {
            preset = minimap2_preset()
            """
            # head closes the pipe once it has enough reads, which gives zcat SIGPIPE (141) under pipefail
            zcat -f ${fastq} | head -n ${params.host_eval_reads * 4} > small.fq || [ \$? -eq 141 ]
            # the reads which are mapped, so that the others are not counted as unmapped
            awk 'NR%4==1 {print substr(\$1,2)}' small.fq > "${unique_id}.${method}.host_sampled.txt"
            minimap2 -ax ${preset} ${refs} small.fq --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > "${unique_id}.${method}.host.sam"
            """
        } else {
            """
            touch "${unique_id}.${method}.host.sam" "${unique_id}.${method}.host_sampled.txt"
            """
        }
}

process extract_microbial_host_hits {

    label "process_medium"
    tag "${unique_id}"
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"

    input:
    tuple val(unique_id), val(method), path(sam_file)
    path ref_bed

    output:
    tuple val(unique_id), val(method), path("query.fasta")

    script:
    """
    samtools view -S -b ${sam_file} | samtools sort - -o ${sam_file.baseName}.sorted.bam
    samtools index ${sam_file.baseName}.sorted.bam
    samtools view -L ${ref_bed} ${sam_file.baseName}.sorted.bam -b -o out.bam
    samtools fasta out.bam > "query.fasta"
    """
}

process blastn_microbial_host_hits {

    label "process_medium_plus_mem"
    tag "${unique_id}"
    label "process_long"
    conda "bioconda::blast=2.16.0"
    container "ncbi/blast"

    input:
    tuple val(unique_id), val(method), path(fasta_file)
    path(blast_db)

    output:
    tuple val(unique_id), val(method), path("results_blastn.txt")

    script:
    if (params.blast_db){
        db = "${blast_db}/nt -num_threads 4"
    } else {
        db = "nt -remote"
    }
    // Best hit is chosen by bitscore downstream, so keep enough targets that human hits do not crowd out the rest.
    // The second search excludes human (taxid 9606) to find the best non-human hit; every line is prefixed
    // with the search it came from ("all" or "nonhuman").
    nonhuman = ""
    if (params.blast_nonhuman) {
        nonhuman = """
      blastn -query ${fasta_file} \\
        -db ${db} \\
        -out nonhuman.txt \\
        -negative_taxids 9606 \\
        -evalue 1e-6 \\
        -perc_identity 85 \\
        -max_target_seqs ${params.blast_max_targets} \\
        -outfmt "6 qseqid sacc sscinames staxids sstart send evalue pident length bitscore qcovs"
      awk 'BEGIN{FS=OFS="\\t"} {print "nonhuman", \$0}' nonhuman.txt >> results_blastn.txt
        """
    }
    """
    : > results_blastn.txt
    if [ -s ${fasta_file} ]; then
      blastn -query ${fasta_file} \\
        -db ${db} \\
        -out all.txt \\
        -evalue 1e-6 \\
        -perc_identity 85 \\
        -max_target_seqs ${params.blast_max_targets} \\
        -outfmt "6 qseqid sacc sscinames staxids sstart send evalue pident length bitscore qcovs"
      awk 'BEGIN{FS=OFS="\\t"} {print "all", \$0}' all.txt >> results_blastn.txt
      ${nonhuman}
    fi
    """
}

process cat_all_microbial_sam_files {

    label "process_low"
    tag "${unique_id}"
    container "community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a"
    publishDir "${params.outdir}/intermediate/", mode: 'copy', pattern: "*.sam"

    input:
    tuple val(unique_id), val(method), path(sam_files)

    output:
    tuple val(unique_id), path("${unique_id}.${method}.microbial.sam")

    script:
    """
    cat \$(ls *.mmp.sam | head -n1) > "${unique_id}.${method}.microbial.sam"
    for sam in \$(ls *.mmp.sam | tail -n+2)
      do
        cat \$sam | tail -n+34 >> "${unique_id}.${method}.microbial.sam"
      done
    """
}


workflow verify_microbial_host_hits {
    take:
        microbial_fastq_ch

    main:
        blast_ch = Channel.empty()
        ref_bed = file("$projectDir/${params.ref_bed}", type: "file", checkIfExists:true)
        extract_microbial_host_hits(microbial_fastq_ch, ref_bed)
        
        if (params.blast_db)
            blast_db = file(params.blast_db, type: "dir", checkIfExists:true)
        else
            blast_db = file("$projectDir/${params.ref_bed}", type: "file", checkIfExists:true) // any file will do to not block

        blastn_microbial_host_hits(extract_microbial_host_hits.out, blast_db)
        blastn_microbial_host_hits.out.collectFile (storeDir: "${params.outdir}/intermediate") { unique_id, method, result -> ["${unique_id}.${method}.results_blastn.txt", result.text]}
                                      .map { f -> [f.simpleName, f] }
                                      .set{ blast_ch }

    emit:
        blast_ch

}

process evaluate_summary {

    label "process_low"
    tag "${unique_id}"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    publishDir "${params.outdir}/${unique_id}/", mode: 'copy'

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
#!/usr/bin/env nextflow

process bam_to_fastq {
    label "process_medium"
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"
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

process minimap2_microbial_chunk {

    label "process_medium"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), val(method), path("microbial.${fastq.baseName}.mmp.sam")
    script:
        if ( params.evaluate_microbial ) {
            preset = ""
            if ( params.read_type == "illumina") {
                preset = "sr"
            } else {
                preset = "map-ont"
            }
            """
            minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > microbial.${fastq.baseName}.mmp.sam
            """
        } else {
            """
            touch "microbial.${method}.mmp.sam"
            """
        }
        
}

workflow minimap2_microbial {
    take:
    fastq_ch
    refs

    main:
    fastq_ch.map{unique_id, method, fastq -> [unique_id, method, fastq.splitFastq(by: params.chunk_size, file:true)]}
            .transpose()
            .set{chunked_fastq_ch}
    minimap2_microbial_chunk(chunked_fastq_ch, refs)
    minimap2_microbial_chunk.out.set{ chunk_sam_ch }
    minimap2_microbial_chunk.out.groupTuple(by: [0,1]).set{collected_sam}
    cat_all_microbial_sam_files(collected_sam)
    cat_all_microbial_sam_files.out.set{ sam_ch }
    
    emit:
    chunk_sam_ch
    sam_ch
}

process minimap2_host {

    label "process_medium"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    //publishDir "${params.outdir}/intermediate/", mode: 'copy', pattern: "*.sam"

    input:
        tuple val(unique_id), val(method), path(fastq)
        path refs
    output:
        tuple val(unique_id), path("${unique_id}.${method}.host.sam")
    script:
        if ( params.evaluate_host == true ) {
            preset = ""
            if ( params.read_type == "illumina") {
                preset = "sr"
            } else {
                preset = "map-ont"
            }
            """
            head -n1000000 ${fastq} > small.fq
            minimap2 -ax ${preset} ${refs} small.fq --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > "${unique_id}.${method}.host.sam"
            """
        } else {
            """
            touch "${unique_id}.${method}.host.sam"
            """
        }
}

process extract_microbial_host_hits {

    label "process_medium"
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"

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
    """
    if [ -s ${fasta_file} ]; then
      blastn -query ${fasta_file} \
        -db ${db} \
        -out results_blastn.txt \
        -evalue 1e-6 \
        -perc_identity 85 \
        -max_target_seqs 5 \
        -outfmt "6 qseqid sacc sscinames staxids sstart send evalue pident length"
    else
      touch "results_blastn.txt"
    fi
    """
}

process cat_all_microbial_sam_files {

    label "process_low"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"
    //publishDir "${params.outdir}/intermediate/", mode: 'copy', pattern: "*.sam"

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
        blastn_microbial_host_hits.out.collectFile { unique_id, method, result -> ["${unique_id}.${method}.results_blastn.txt", result.text]}
                                      .map { f -> [f.simpleName, f] }
                                      .set{ blast_ch }

    emit:
        blast_ch

}

process evaluate_summary {

    label "process_low"
    container 'community.wave.seqera.io/library/simplesam_numpy_pandas_pip_taxoniq:3af1649bdaf86fca'
    publishDir "${params.outdir}/${unique_id}/", mode: 'copy'

    input:
    tuple val(unique_id), val(classifier), path(report), path(host_sam), path(microbial_sam), path(blast_result)

    output:
    path "${unique_id}*_summary.csv", emit: summary
    path "${unique_id}*_full.csv", emit: full
    path "${unique_id}*_data.csv", emit: data
    path "${unique_id}*_taxa.csv", emit: taxa, optional:true
    path "${unique_id}*_accs.csv", emit: accs, optional:true

    script:
    """
    evaluate.py \
      -i ${report} \
      --microbial_sam ${microbial_sam} \
      --host_sam ${host_sam} \
      --blast_result ${blast_result} \
      -p "${unique_id}_${classifier}"
    """
}
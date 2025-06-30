#!/usr/bin/env nextflow

process bam_to_fastq {
    label "process_medium"
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"

    input:
    tuple val(unique_id), path(bam)

    output:
    tuple val(unique_id), path("${bam.baseName}.fastq")

    script:
    """
    samtools fastq ${bam} > "${bam.baseName}.fastq"
    """
}

process minimap2_microbial {

    label "process_medium"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(fastq)
        path refs
    output:
        tuple val(unique_id), path("microbial.mmp.sam")
    script:
        if ( params.evaluate_microbial ) {
            preset = ""
            if ( params.read_type == "illumina") {
                preset = "sr"
            } else {
                preset = "map-ont"
            }
            """
            minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > microbial.mmp.sam
            """
        } else {
            """
            touch "microbial.mmp.sam"
            """
        }
        
}

process minimap2_host {

    label "process_medium"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(fastq)
        path refs
    output:
        tuple val(unique_id), path("host.mmp.sam")
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
            minimap2 -ax ${preset} ${refs} small.fq --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > host.mmp.sam
            """
        } else {
            """
            touch "host.mmp.sam"
            """
        }
}

process extract_microbial_host_hits {

    label "process_medium"
    conda "bioconda::samtools=1.21"
    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"

    input:
    tuple val(unique_id), path(sam_file)
    path ref_bed

    output:
    tuple val(unique_id), path("query.fasta")

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
    conda "bioconda::blast=2.16.0"
    container "ncbi/blast"

    input:
    tuple val(unique_id), path(fasta_file)
    path(blast_db)

    output:
    tuple val(unique_id), path("results_blastn.txt")

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
        -perc_identity 80 \
        -max_target_seqs 5 \
        -outfmt "6 qseqid sacc sscinames staxids sstart send evalue pident length"
    else
      touch "results_blastn.txt"
    fi
    """
}

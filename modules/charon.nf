#!/usr/bin/env nextflow

process run_charon {

    label "process_medium_plus_mem"
    container 'docker.io/rmcolq/charon:v1.0.5'

    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("charon_${unique_id}_microbial.fq.gz"), emit: microbial_fastq
    tuple val(unique_id), path("charon_${unique_id}_human.fq.gz"), emit: human_fastq
    tuple val(unique_id), path("charon_${unique_id}.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      --log charon_${unique_id}.log \
      --extract all \
      --prefix charon_${unique_id} \
      -t ${task.cpus} \
      --min_length 20 > charon_${unique_id}.out
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
            minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > host.mmp.sam
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

    label "process_medium"
    conda "bioconda::blast=2.16.0"
    container "ncbi/blast"

    input:
    tuple val(unique_id), path(fasta_file)

    output:
    tuple val(unique_id), path("results_blastn.txt")

    script:
    """
    blastn -query ${fasta_file} \
      -db ${params.blast_db} \
      -out results_blastn.txt \
      -evalue 1e-6 \
      -perc_identity 90 \
      -outfmt "6 qseqid sacc sscinames staxids evalue pident length"
    """
}

process evaluate_summary {

    container 'community.wave.seqera.io/library/simplesam_numpy_pandas:ea9b7172ad7bff36'
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(unique_id), path(charon_report), path(host_sam), path(microbial_sam), path(blast_result)

    output:
    path "${unique_id}_summary.csv", emit: summary
    path "${unique_id}*_data.csv", emit: data

    script:
    """
    evaluate.py \
      -i ${charon_report} \
      --microbial_sam ${microbial_sam} \
      --host_sam ${host_sam} \
      --blast_result ${blast_result} \
      -p "${unique_id}"
    """
}


workflow evaluate_charon {
    take:
        fastq_ch
    main:

    db = file(params.db, type: "file", checkIfExists:true)
    refs = file("$projectDir/${params.refs}", type: "file", checkIfExists:true)

    run_charon(fastq_ch, db)
    minimap2_microbial(run_charon.out.microbial_fastq, refs)
    minimap2_host(run_charon.out.human_fastq, refs)

    if ( params.evaluate_microbial ){
        blast_ch = Channel.empty()
        ref_bed = file("$projectDir/${params.ref_bed}", type: "file", checkIfExists:true)
        extract_microbial_host_hits(minimap2_microbial.out, ref_bed)
        blastn_microbial_host_hits(extract_microbial_host_hits.out)
        blastn_microbial_host_hits.out.set{ blast_ch }
    } else {
        blast_ch = Channel.empty()
    }

    run_charon.out.result
             .combine(minimap2_host.out, by: 0)
             .combine(minimap2_microbial.out, by: 0)
             .combine(blast_ch, by: 0)
             .set{ eval_ch }

    evaluate_summary(eval_ch)

}

#!/usr/bin/env nextflow

process run_charon_microbial {

    label "process_high"
    container 'docker.io/rmcolq/charon:latest'

    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("charon_${unique_id}_microbial.fq.gz"), emit: fastq
    tuple val(unique_id), path("charon_${unique_id}_microbial.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      --log charon_${unique_id}_microbial.log \
      --extract microbial \
      --extract_file charon_${unique_id}_microbial.fq.gz \
      -t ${task.cpus} \
      --min_length 20 > charon_${unique_id}_microbial.out
    """
}

process run_charon_host {

    label "process_high"
    container 'docker.io/rmcolq/charon:latest'

    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("charon_${unique_id}_host.fq.gz"), emit: fastq
    tuple val(unique_id), path("charon_${unique_id}_host.out"),  emit: result

    script:
    """
    charon dehost ${fastq} \
      --db ${db} \
      --log charon_${unique_id}_host.log \
      --extract human \
      --extract_file charon_${unique_id}_host.fq.gz \
      -t ${task.cpus} \
      --min_length 20 > charon_${unique_id}_host.out
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
        preset = ""
        if ( params.read_type == "illumina") {
            preset = "sr"
        } else {
            preset = "map-ont"
        }
        """
        minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > microbial.mmp.sam
        """
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
        preset = ""
        if ( params.read_type == "illumina") {
            preset = "sr"
        } else {
            preset = "map-ont"
        }
        """
        minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > host.mmp.sam
        """
}

process evaluate_summary {

    container 'community.wave.seqera.io/library/simplesam_numpy_pandas:ea9b7172ad7bff36'

    input:
    tuple val(unique_id), path(charon_report), path(host_sam), path(microbial_sam)

    output:
    path "${unique_id}_summary.csv", emit: summary
    path "${unique_id}_data.csv", emit: data

    publishDir "${params.outdir}/", mode: 'copy'

    script:
    """
    evaluate.py \
      -i ${charon_report} \
      --microbial_sam ${microbial_sam} \
      --host_sam ${host_sam} \
      -p "${unique_id}"
    """
}


workflow evaluate_charon {
    unique_id = "${params.unique_id}"
    fastq = file(params.fastq, type: "file", checkIfExists:true)
    //report = file(params.charon_report, type: "file", checkIfExists:true)
    fastq_ch = Channel.from([[unique_id, fastq]])
    //report_ch = Channel.from([[unique_id, report]])

    db = file(params.db, type: "file", checkIfExists:true)
    refs = file("$projectDir/resources/T2T_ref_plus_EBV.fa.gz", type: "file", checkIfExists:true)

    run_charon_microbial(fastq_ch, db)
    minimap2_microbial(run_charon_microbial.out.fastq, refs)

    run_charon_host(fastq_ch, db)
    minimap2_host(run_charon_host.out.fastq, refs)

    run_charon_microbial.out.result
             .combine(minimap2_host.out, by: 0)
             .combine(minimap2_microbial.out, by: 0)
             .set{ eval_ch }

    evaluate_summary(eval_ch)

}

#!/usr/bin/env nextflow

process run_charon_microbial {

    label "process_medium"
    container 'community.wave.seqera.io/library/pip_numpy_pandas:426ad974eac1c1db'

    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("charon_${unique_id}_microbial.fq.gz"), emit: fastq
    tuple val(unique_id), path("charon_${unique_id}_microbial.log"), path("charon_${unique_id}_microbial.out"),  emit: result

    script:
    """
    charon classify ${fastq} \
      --db ${db} \
      --log charon_${unique_id}_microbial.log \
      --extract microbial \
      --extract_file charon_${unique_id}_microbial.fq.gz \
      -t ${task.cpu} \
      --min_length 20 > charon_${unique_id}_microbial.out
    """
}

process run_charon_human {

    label "process_medium"
    container 'community.wave.seqera.io/library/pip_numpy_pandas:426ad974eac1c1db'

    input:
    tuple val(unique_id), path(fastq)
    path(db)

    output:
    tuple val(unique_id), path("charon_${unique_id}_microbial.fq.gz"), emit: fastq
    tuple val(unique_id), path("charon_${unique_id}_microbial.log"), path("charon_${unique_id}_microbial.out"),  emit: result

    script:
    """
    charon classify ${fastq} \
      --db ${db} \
      --log charon_${unique_id}_human.log \
      --extract human \
      --extract_file charon_${unique_id}_human.fq.gz \
      -t ${task.cpu} \
      --min_length 20 > charon_${unique_id}_human.out
    """
}

process minimap2_microbial {

    label "process_low"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(fastq)
        path refs
    output:
        tuple val(unique_id), val(database_name), path(kreport), path(reads), path("microbial.mmp.sam")
    script:
        preset = ""
        if ( params.read_type == "illumina") {
            preset = "sr"
        } else {
            preset = "map-ont"
        }
        """
        minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > hcid.mmp.sam
        """
}

process minimap2_human {

    label "process_low"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(fastq)
        path refs
    output:
        tuple val(unique_id), val(database_name), path(kreport), path(reads), path("human.mmp.sam")
    script:
        preset = ""
        if ( params.read_type == "illumina") {
            preset = "sr"
        } else {
            preset = "map-ont"
        }
        """
        minimap2 -ax ${preset} ${refs} ${fastq} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > hcid.mmp.sam
        """
}

process evaluate {

    container 'community.wave.seqera.io/library/numpy_pandas_python:707bf65260f82405'

    input:
    tuple val(unique_id), path(charon_report), path(human_sam), path(microbial_sam)

    output:
    path "${unique_id}_summary.csv"

    publishDir "${params.outdir}/", mode: 'copy'

    script:
    """
    evaluate.py \
      -i ${charon_report} \
      --microbial_sam ${microbial_sam} \
      --host_sam ${host_sam} \
      --output "${unique_id}_summary.csv"
    """
}


workflow evaluate_charon {
    unique_id = "${params.unique_id}"
    fastq = file(params.fastq, type: "file", checkIfExists:true)
    report = file(params.charon_report, type: "file", checkIfExists:true)
    fastq_ch = Channel.from([[unique_id, fastq]])
    report_ch = Channel.from([[unique_id, report]])

    db = file(params.db, type: "file", checkIfExists:true)
    refs = file("$projectDir/resources/hcid_refs.fa.gz", type: "file", checkIfExists:true)

    run_charon_microbial(fastq_ch, db)
    minimap2_microbial(run_charon_microbial.out.fastq, refs)

    run_charon_human(fastq_ch, db)
    minimap2_human(run_charon_human.out.fastq, refs)

    report_ch.combine(minimap2_microbial.out, by: 0)
             .combine(minimap2_human.out, by: 0)
             .set{ eval_ch }

    evaluate(eval_ch)

}

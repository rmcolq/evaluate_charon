#!/usr/bin/env nextflow

//take input as --reports and--metadata

/*
 * A Python script which parses the output of the previous script
 */
process run_charon {

    label "process_long"
    container 'community.wave.seqera.io/library/pip_numpy_pandas:426ad974eac1c1db'

    input:
    tuple val(unique_id), path(fastq)

    output:
    path "text_files"

    script:
    """
    make_r_files.py --reports ${reports.join(' ')} --metadata ${metadata} --site_key ${site_key} --output_dir text_files/
    """
}


/*
 * An R script which produces output for shannon's diversity graphs
 */

process get_shannon_plot {

    container 'community.wave.seqera.io/library/r-argparse_r-crayon_r-dplyr_r-ggplot2_r-vegan:eb552a73894bf74c'
    input:
    path text_files
    output:
    path "plots"

    script:
    """
    make_r_plots.R ${text_files}/ plots/
    """
}

/*
 * A python script which produces output for making the html report
 */

process make_report {

    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:44e99f335376fa3b'

    input:
    path reports
    path metadata
    path site_key
    path template
    path plots

    output:
    path "summary_report/*.html"

    publishDir "${params.outdir}/", mode: 'copy' // Publish final report to local directory specified in params.config

    script:
    """
    make_sum_report.py \
      --reports ${reports.join(' ')} \
      --metadata ${metadata} \
      --site_key ${site_key} \
      --plots_dir ${plots}/ \
      --final_report summary_report/ \
      --template ${template}
    """
}


workflow evaluate_negative_controls {
    report_list = params.reports?.split('\n') as List
    Channel
        .fromPath(report_list)
        .flatten()
        .collect()
        .set { reports }
    reports.view()


    metadata_file = file(params.metadata, type: "file", checkIfExists:true)
    Channel
        .fromPath(metadata_file)
        .set { metadata }

    site_key = file(params.site_key, type: "file", checkIfExists:true)

    template = file("$baseDir/bin/summary_report_template.html")

    make_shannon_script(reports, metadata, site_key)
    get_shannon_plot(make_shannon_script.out)
    make_report(reports, metadata, site_key, template, get_shannon_plot.out)
    println "Report will be generated in ${params.outdir}"
}

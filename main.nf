include { evaluate_negative_controls } from './modules/summary_report'

// Check if required parameters are provided
if (!params.reports) {
    exit 1, "Please provide --reports when running Nextflow."
}

if (!params.metadata) {
    exit 1, "Please provide --metadata when running Nextflow."
}

workflow {
    evaluate_negative_controls()
}
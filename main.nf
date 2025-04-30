include { evaluate_charon } from './modules/charon'

// Check if required parameters are provided
if (!params.unique_id) {
    exit 1, "Please provide --unique_id when running Nextflow."
}

if (!params.fastq) {
    exit 1, "Please provide --fastq when running Nextflow."
}

if (!params.db) {
    exit 1, "Please provide --db when running Nextflow."
}

if (params.read_type == "illumina") {
    exit 1, "I don't support paired end reads yet"
}

workflow {
    evaluate_charon()
}
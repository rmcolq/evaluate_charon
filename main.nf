include { evaluate_charon } from './modules/charon'



workflow {

    // Check if required parameters are provided
    if (!params.db) {
        exit 1, "Please provide --db when running Nextflow."
    }

    if (params.read_type == "illumina") {
        exit 1, "I don't support paired end reads yet"
    }

    if (params.fastq_dir){
        fastq_ch = Channel.fromPath("${params.fastq_dir}/*.f*q*", type: "file", checkIfExists: true).map { [it.baseName.replace(".fq","").replace(".fastq",""), it] }
        fastq_ch.view()
        evaluate_charon(fastq_ch)
    } else if (params.fastq){
        unique_id = "${params.unique_id}"
        fastq = file(params.fastq, type: "file", checkIfExists:true)
        fastq_ch = Channel.from([[unique_id, fastq]])
        evaluate_charon(fastq_ch)
    } else {
        exit 1, "Please either provide a directory of fastq with --fastq_dir (each file for sample) or --fastq and --unique_id"
    }

}
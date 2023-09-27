#!/usr/bin/env nextflow

/// Specify to use Nextflow DSL version 2
nextflow.enable.dsl=2

/// Import modules and subworkflows
include { preprocessing } from './subworkflows/local/preprocessing.nf'
include { quantification } from './subworkflows/local/quantification.nf'

// Log the parameters
log.info """\

=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||                        INSERT PIPELINE NAME                             
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  Parameters                                                             
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
||  input_dir   : ${params.input_dir}                                     
||  outDir      : ${params.output_dir}                                        
||  workDir     : ${workflow.workDir}                                     
=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

"""
// Help Message to prompt users to specify required parameters
def help() {
    log.info"""
  Usage:  nextflow run main.nf --input_dir <path_to_fastq_dir> 

  Required Arguments:

  --input    Path to directory containing fastq files.

  Optional Arguments:

  --outDir	Path to output directory. 
	
""".stripIndent()
}

/// Define the main workflow
workflow {
    // println "$params.samplesheet"
    // fastq_ch = Channel.fromPath("${params.samplesheet}").splitCsv(header: true)
    // fastq_ch.view()
    // /// Define the input channels
    // fastq_ch = Channel.fromFilePairs("${params.input_dir}/*{1,2}.fastq.gz", size: 2)
    //                    .ifEmpty { exit 1, "No fastq files found in ${params.input_dir}" }
    
    fastq_ch = Channel.fromPath("${params.input_dir}/*.fastq.gz")
                        .ifEmpty { exit 1, "No fastq files found in ${params.input_dir}" }
    fastq_ch.view()

    /// Run the subworkflow
    preprocessed_files = preprocessing(fastq_ch)
    preprocessed_files.view()
    //quantification(preprocessed_files)  // How do I specify the output?
}

workflow.onComplete {
    log.info "Pipeline completed at: ${new Date().format('dd-MM-yyyy HH:mm:ss')}"
}

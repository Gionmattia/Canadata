

process FASTQC {

    label 'cpu_12'
	publishDir "${params.output_dir}/fastqc", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
    // FASTQC can work with two files, but the output generation is wrong:
    // It tries to save both files together.
    // Solution: process the files singularly at this stage.


	input:
        path(fastq) // tuple val(sample_id), file(fastq)

	output:
	    path "*_fastqc.html", emit: fastqc_html

        tuple path(fastq), path("${fastq.baseName.minus('.fastq')}_fastqc/fastqc_data.txt"), emit: fastqc_data
        //path "${fastq.baseName.minus('.fastq')}_fastqc/fastqc_data.txt", emit: fastqc_data_txt


    script:
        """
        fastqc --extract -q $fastq --adapters ${projectDir}/scripts/adapter_list.tsv --dir ${projectDir}/data/temp
        """
}
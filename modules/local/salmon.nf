

process SALMON_QUANT {

    tag 'medium'

	publishDir "${params.study_dir}/salmon", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
	    file fastq 
        file index

	output:
	    path "*_fastqc.html", emit: fastqc_html
        path "${fastq.baseName}_fastqc/fastqc_data.txt", emit: fastqc_data

    script:
        """
        salmon quant \
            -i $index \
            -p 10 \
            -l A \
            -1 ${input_dir}/"$sample""_1_less_rRNA""$extension" \
            -2 ${input_dir}/"$sample""_2_less_rRNA""$extension" \
            -o "$sample""quant"        
    """
}
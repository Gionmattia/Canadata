

process SALMON_QUANT {

    tag 'medium'
	publishDir "${params.output_dir}/salmon_quant", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
	    tuple val(sample_id), path("*.fastq.gz") from no_rRNA_fastq_ch // OR COULD BE "from less_rRNA.no_rRNA_fastq_ch"
        // 
	output:

        path "${fastq.baseName}_fastqc/fastqc_data.txt", emit: fastqc_data

    // Need to check the -1 and -2 option actually refers the naming convetions applied so far. It could be different from
    // what expected

    script:
        """
        salmon quant \
            -i ${params.salmon_index} \
            -p 10 \
            -l A \
            -1 ${input_dir}/"$sample""_1_less_rRNA""$extension" \
            -2 ${input_dir}/"$sample""_2_less_rRNA""$extension" \
            -o "$sample""quant"        
    """
}

// -o should make use of "sampleID" to name the stuff.
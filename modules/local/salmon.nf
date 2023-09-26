

process SALMON_QUANT {

    // THESE INFO SHOULD ALL BE PUT WITHIN THE nextflow.config FILE (or so I suppose...it just makes sense)
    tag 'medium'
    container "${params.singularity_path}"  // Path to your local Singularity image file
    containerOptions "--bind ./data"  // Optional: bind specific paths if needed

    publishDir "${params.output_dir}/salmon_quant", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
	    tuple val(sample_id), path("*.fastq.gz")// from no_rRNA_fastq // OR COULD BE "from less_rRNA.no_rRNA_fastq_ch"

	output:

        path "${fastq.baseName}_salmon_quants.sf", emit: salmon_counts

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
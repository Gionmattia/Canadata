process TXIMPORT {

    tag 'medium'
    publishDir "${params.output_dir}/gene_counts", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
        path (quant_file)

	output:
    
        path ("*.tsv")

    script:
    
    """
    Rscript ${projectDir}/scripts/tximport_convert.R ${quant_file} ${params.annotation_path} 
    """

}

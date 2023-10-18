process TXI_CONVERSION {

    tag 'medium'
    publishDir "${params.output_dir}/txi_gene_counts", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
        path (quant_file)

	output:
    
    //    path ("*/*_quant.sf"), emit: quant_sf define after script block

    script:
    
    """
    rstudio-something 
    """


}

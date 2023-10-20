process SALMON_QUANT {
    label 'cpu_12'
    //tag 'medium'
    publishDir "${params.output_dir}/salmon_quants", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
        tuple path(reads1), path(reads2)

	output:
        path ("*/*")
        path ("*/*_quant.sf"), emit: quant_sf

    script:

    def sample_path = "${reads1}"
    def sample_name = sample_path.split("/")[-1]
    def sample_id = "${sample_name.minus("_less_rRNA_1.fastq")}"

    """
    singularity exec -B ${projectDir} ${projectDir}/singularity/salmon:1.10.1--h7e5ed60_0 salmon quant -i ${params.salmon_index} \
    -p 10 \
    -l A \
    -1 ${reads1} \
    -2 ${reads2} \
    -o "${sample_id}"

    mv "${sample_id}/quant.sf" "${sample_id}/${sample_id}_quant.sf"
    """


}

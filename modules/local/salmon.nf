process SALMON_QUANT {

    tag 'medium'
    publishDir "${params.output_dir}/salmon_quants", mode: 'copy'

    //container "/data2/Canadata/Canadata/singularity/salmon:1.10.1--h7e5ed60_0"
    //runOptions "--bind ${projectDir}"

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
	input:
	    tuple val(sample_id), path("*.fastq.gz")// from no_rRNA_fastq // OR COULD BE "from less_rRNA.no_rRNA_fastq_ch"

	output:

        //path "${params.output_dir}/salmon_quants/${sample_id}_salmon_quants.sf", emit: salmon_counts

    // Issues:
    // 1 The code does not allow for relative paths. Absolute only.
    // 2) The FASTQ files are asynchronous. Need to find a better adpater remover or use the options more sensibly.

    script:
    """
    singularity exec -B /data2/Canadata/Canadata --workdir /data2/Canadata/Canadata /data2/Canadata/Canadata/singularity/salmon:1.10.1--h7e5ed60_0 salmon quant -i ${params.salmon_index} \
    -p 10 \
    -l A \
    -1 data/output/bowtie/"${sample_id}__1.fastq_clipped_less_rRNA.fastq.gz" \
    -2 /data2/Canadata/Canadata/data/output/bowtie/"${sample_id}__2.fastq_clipped_less_rRNA.fastq.gz" \
    -o "${sample_id}"
    """
}

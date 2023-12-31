process FIND_ADAPTERS {
    label 'cpu_12'
	publishDir "${params.output_dir}/adapter_reports", mode: 'copy'
    
    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        tuple path(raw_fastq), path(fastqc_data_txt) // tuple val(sample_id), file(raw_fastq)
        //path(fastqc_data)

    output:
        file "*_adapter_report.fa"
        // The following is a dummy string used to create a dependency between this process and FASTP
        val "synchronising node", emit: dummy_string

    script:
        """
        python3 ${projectDir}/scripts/get_adapters.py -i $fastqc_data_txt -a ${projectDir}/scripts/adapter_list.tsv -o "${raw_fastq}_adapter_report.fa"
        echo "synchronising node"
        """
}
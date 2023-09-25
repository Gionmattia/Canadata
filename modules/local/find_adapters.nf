process FIND_ADAPTERS {
	publishDir "${params.study_dir}/adapter_reports", mode: 'copy'
    
    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }

    input:
        tuple val(sample_id), file(raw_fastq)
        file fastqc_data

    output:
        file "*_adpater_report.fa"

    script:
        """
        python3 $projectDir/scripts/get_adapters.py -i $fastqc_data -a $projectDir/scripts/adapter_list.tsv -o "${raw_fastq}_adpater_report.fa"
        """
}
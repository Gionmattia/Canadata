process FASTP {
    publishDir "${params.output_dir}/fastp", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/fastp", mode: 'copy', pattern: '*.html'

    // errorStrategy 'ignore'
    

    // Currently takes as input the tuple, with the paired files.
    input:
    tuple val(sample_id), file(raw_fastq)  // path(raw_fastq) 

    
    // Change the output so it gets the two files in output
    output:
    tuple path ("*${sample_id}*1.fastq_clipped.fastq"), path ("*${sample_id}*2.fastq_clipped.fastq"), emit: trimmed_fastq

    //FASTP produces only one report, even if there are two read-files
    path ("*${sample_id}*.json"), emit: fastp_json
    path ("*${sample_id}*.html"), emit: fastp_html


 script:

    def reads1 = raw_fastq[0]
    def reads2 = raw_fastq[1]

    def file1 = "${projectDir}/data/output/adapter_reports/${sample_id}__1.fastq.gz_adpater_report.fa"
    def file2 = "${projectDir}/data/output/adapter_reports/${sample_id}__2.fastq.gz_adpater_report.fa"
    
    """
    bash ${projectDir}/scripts/FASTP_run.sh "$reads1" "$reads2" "$file1" "$file2"
    """
}

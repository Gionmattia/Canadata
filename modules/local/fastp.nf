process FASTP {
    publishDir "${params.output_dir}/fastp", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/fastp", mode: 'copy', pattern: '*.html'

    // errorStrategy 'ignore'
    

    // Currently takes as input the tuple, with the paired files.
    input:
    tuple val(sample_id), file(raw_fastq)  // path(raw_fastq) 

    // Need to take in BOTH adapter reports!
    // OR
    // Could just take the reports from the folder where they are stored instead, so that I don't
    // have the fuss of moving them across.
    
    path(adapter_report)  // file adapter_report

    output:
    path ("${raw_fastq.baseName}_clipped.fastq"), emit: trimmed_fastq
    path '*.json', emit: fastp_json
    path '*.html', emit: fastp_html

	
	script: 
	"""
    fastp \
    -i $raw_fastq \
    -o ${raw_fastq.baseName}_clipped.fastq \
    --length_required 20 \
    --adapter_fasta $adapter_report \
    --json ${raw_fastq.baseName}.json \
    --html ${raw_fastq.baseName}.html \
    """
}

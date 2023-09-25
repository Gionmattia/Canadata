process FASTP {
    publishDir "$params.study_dir/fastp", mode: 'copy', pattern: '*.json'
    publishDir "$params.study_dir/fastp", mode: 'copy', pattern: '*.html'

    // errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), file(raw_fastq)
    file adapter_report

    output:
    path '*_clipped.fastq', emit: trimmed_fastq
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

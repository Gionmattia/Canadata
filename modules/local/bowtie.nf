process rRNA_REMOVAL {

    // This processes removes rRNAs, Mt_rRNAs and Mt_tRNAs

    label 'cpu_12'  // custom label for requesting 12 CPUs, define in your resources configuration
    publishDir "${params.output_dir}/bowtie", mode: 'copy' // NOT SURE I NEED THIS LINE OF CODE, STUFF COULD BE EMITTED IN A CHANNEL

    input:
    path(trimmed_fastq)  // tuple val(sample_id), path(trimmed_fastq)

    output:
    // Tried to implement the same code as in the other modules, so the base_name is used.
    path("${trimmed_fastq.baseName}_less_rRNA.fastq.gz") // , emit: no_rRNA_fastq  //originally it was ("${sample_id}_less_rRNA.fastq.gz")

    script:
    """
    bowtie -p ${task.cpus} -v 3 ${params.rRNA_index} -q ${trimmed_fastq} --un ${trimmed_fastq.baseName}_less_rRNA.fastq.gz >> screen_output.txt 2>&1
    """
}

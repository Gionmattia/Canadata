process rRNA_REMOVAL {

    // This processes removes rRNAs, Mt_rRNAs and Mt_tRNAs

    label 'cpu_12'  // custom label for requesting 12 CPUs, define in your resources configuration
    publishDir "${params.output_dir}/bowtie", mode: 'copy'

    input:
    // NEED TO UPDATE TO RECENT VERSION!

    tuple path(sample_id), path(read1), path(read2)
    // tuple path("*1.fastq_clipped.fastq"), path("*2.fastq_clipped.fastq")


    output:
    // Emit as a tuple for salmon later on?
    tuple path("${read1.baseName}_less_rRNA.fastq.gz"), path("${read2.baseName}_less_rRNA.fastq.gz"), emit: no_rRNA_fastq
    //path("${trimmed_fastq.baseName}_less_rRNA.fastq.gz") // , emit: no_rRNA_fastq 

    script:

    // NEED TO IMPLEMENT IN THE CODE THE PE MODE (-1,-2 args?)
    script:
    """
    bowtie -p ${task.cpus} -v 3 ${params.rRNA_index} -1 ${read1} -2 ${read2} --un ${sample_id}_less_rRNA.fastq.gz >> screen_output.txt 2>&1
    """



    //"""
    //bowtie -p ${task.cpus} -v 3 ${params.rRNA_index} -q ${trimmed_fastq} --un ${trimmed_fastq.baseName}_less_rRNA.fastq.gz >> screen_output.txt 2>&1
    //"""
}

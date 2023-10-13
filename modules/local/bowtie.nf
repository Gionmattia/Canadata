process rRNA_REMOVAL {

    // This processes removes rRNAs, Mt_rRNAs and Mt_tRNAs

    label 'cpu_12'  // custom label for requesting 12 CPUs, define in your resources configuration
    publishDir "${params.output_dir}/bowtie", mode: 'copy'

    input:
    // NEED TO UPDATE TO RECENT VERSION!

    tuple path(sample_id), path(reads)


    output:
    // CHECK AND FIX THIS
    tuple path("${read1.baseName}_less_rRNA.fastq.gz"), path("${read2.baseName}_less_rRNA.fastq.gz"), emit: no_rRNA_fastq
    //path("${trimmed_fastq.baseName}_less_rRNA.fastq.gz") // , emit: no_rRNA_fastq 


    script:

    // NEED TO IMPLEMENT IN THE CODE THE PE MODE (-1,-2 args?)
    // NEED TO UNPACK THE TUPLE AS USUAL
    """
    bowtie -p ${task.cpus} -v 3 ${params.rRNA_index} -1 ${read1} -2 ${read2} --un ${sample_id}_less_rRNA.fastq.gz >> screen_output.txt 2>&1
    """

    // TO DO ON MONDAY:
    // - Fix cardinality, the tuple needs to be unpacked before giving it to the script.
    // - Updated script code to work on PE
    // - Check how the "un" option names the read files and if consistent with the names so far.
    // - Test run and troubleshoot. There will be minor issues.
    // - Check output formats.
    // PREPARE TO THE TRANSITION TO SALMON. NOW EVERYTHING SHOULD BE SUBMITTED AS A TUPLE CHANNEL, FINALLY.
    // Nb. remember to fix the salmon.nf to use the absolute paths given by ${projectDir}

    //"""
    //bowtie -p ${task.cpus} -v 3 ${params.rRNA_index} -q ${trimmed_fastq} --un ${trimmed_fastq.baseName}_less_rRNA.fastq.gz >> screen_output.txt 2>&1
    //"""
}

process rRNA_REMOVAL {

    // This processes removes rRNAs, Mt_rRNAs and Mt_tRNAs

    label 'cpu_12' 
    publishDir "${params.output_dir}/bowtie", mode: 'copy'

    input:
    tuple path(reads1), path(reads2)

    output:
    // When run in PE mode, bowtie will generate two files for each read-file, with the same
    // <basename> + "_1" or "_2"
    tuple path("*less_rRNA_1.fastq"), path("*less_rRNA_2.fastq"), emit: no_rRNA_fastq

    script:
    // This script just retrieves the sample name from the path, to use as <basename> in bowtie
    // This is done because bowtie does not allow to keep the names of the original read files
    // when used to filter away RNAs...
    
    def str_sample = "${reads1}"
    def sample_name = str_sample.split("/")[-1]
    def stripped_sample_name = "${sample_name.minus("__1.fastq_clipped.fastq")}"
    
    """
    bowtie -p ${task.cpus} -v 3 ${params.rRNA_index} \
    -1 ${reads1} \
    -2 ${reads2} \
    --un ${stripped_sample_name}_less_rRNA.fastq >> screen_output.txt 2>&1
    """

    // PREPARE TO THE TRANSITION TO SALMON. NOW EVERYTHING SHOULD BE SUBMITTED AS A TUPLE CHANNEL, FINALLY.
    // Nb. remember to fix the salmon.nf to use the absolute paths given by ${projectDir}

}

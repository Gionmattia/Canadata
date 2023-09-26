// workflow for counting reads 

include { SALMON_QUANT } from '../../modules/local/salmon.nf'

workflow quantification {

    take: 
        less_rRNA.no_rRNA_fastq_ch    // needs to be determined from the output of the previous subworkflow.

    main:
        salmon_counts       =   SALMON_QUANT          ( less_rRNA.no_rRNA_fastq_ch )
        
    emit:
        // Needs to be determined based on the output of the SALMON_QUANT process

}

// MAIN IDEA: This subworkflow ONLY gets the count-data as a <salmon_quant>.sf file.
// This because the next steps will be more dependent on the experiment itself
// and because this modularity helps me in the future: If I want to use STAR instead, or ORFik,
// I can just input everything through that sub-workflow instead of this one.
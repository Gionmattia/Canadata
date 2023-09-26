// workflow for counting reads 

include { SALMON_QUANT } from '../../modules/local/salmon.nf'

workflow quantification {

    take: 
        no_rRNA_fastq    // OR Coudl be "less_rRNA.no_rRNA_fastq_ch"   needs to be determined from the output of the previous subworkflow.

    main:
        salmon_counts       =   SALMON_QUANT          ( no_rRNA_fastq )
        
    //emit:
        //LOREM  IPSUM     // Needs to be determined based on the output of the SALMON_QUANT process

}

// MAIN IDEA: This subworkflow ONLY gets the count-data as a <salmon_quant>.sf file and its outputs converter.
// If I want to move to STAR or ORFik, I will just need to add a different subworkflow in here.

// This sub-workflow should work as follows:
// INPUT: fastqc_less_rRNA files (in pairs)
// SALMON_QUANT --> .sf output
// TXIMPORT --> converts the .sf output from transcripts to genes
// TSV_MERGER --> extract the target info (gene counts)
// OUTPUT: a .tsv file with all count data for all genes in the samples.

// next this .tsv file should be opened and processed by deltaTE.
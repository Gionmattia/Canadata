// workflow for preprocessing data 

include { FASTQC } from '../../modules/local/fastqc.nf'
include { FIND_ADAPTERS } from '../../modules/local/find_adapters.nf'
include { FASTP } from '../../modules/local/fastp.nf'
include { rRNA_REMOVAL } from '../../modules/local/bowtie.nf'

workflow preprocessing {

    take: 
        //fastq_ch  ++ see below the operation
        fastq_ch_paired

    main:
        // Creates a single-file channel, from the one of pairs. 
       fastq_ch_single = fastq_ch_paired.flatMap { sample_id, files -> files.collect { [it] }}
        
        // NEED TO TRY IMPLEMENT THE fastq_ch_single below!
        fastqc_ch           =   FASTQC          ( fastq_ch_single )
        adapter_ch          =   FIND_ADAPTERS   ( fastq_ch_single, fastqc_ch.fastqc_data )

        // Need to add to its arguments also the second adapter report channel...?
        trimmed_fastq_ch    =   FASTP           ( fastq_ch_paired)  //  ( fastq_ch_paired, adapter_ch )  
        // ORIGINALLY trimmed_fastq_ch    =   FASTP           ( fastq_ch, adapter_ch )
        //less_rRNA           =   rRNA_REMOVAL    ( trimmed_fastq_ch.trimmed_fastq )
        //less_rRNA.collect()

}
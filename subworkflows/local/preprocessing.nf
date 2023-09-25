// workflow for preprocessing data 

include { FIND_ADAPTERS } from '../../modules/local/find_adapters.nf'
include { FASTP } from '../../modules/local/fastp.nf'
include { FASTQC } from '../../modules/local/fastqc.nf'
include { rRNA_REMOVAL } from '../../modules/local/bowtie.nf'

workflow preprocessing {

    take: 
        fastq_ch    

    main:
        fastqc_ch           =   FASTQC          ( fastq_ch )
        adapter_ch          =   FIND_ADAPTERS   ( fastq_ch, fastqc_ch.fastqc_data )
        trimmed_fastq_ch    =   FASTP           ( fastq_ch, adapter_ch )
        // Here starts the new parts of code I implemented
        less_rRNA           =   rRNA_REMOVAL    ( trimmed_fastq_ch )
        
    emit:
        less_rRNA.no_rRNA_fastq_ch   // this should be enough to emit the channel with the sample_ID and the files.
        //trimmed_fastq_ch.trimmed_fastq
}

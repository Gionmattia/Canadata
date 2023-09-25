// workflow for preprocessing data 

include { FIND_ADAPTERS } from '../../modules/local/find_adapters.nf'
include { FASTP } from '../../modules/local/fastp.nf'
include { FASTQC } from '../../modules/local/fastqc.nf'

workflow preprocessing {

    take: 
        fastq_ch    

    main:
        fastqc_ch           =   FASTQC          ( fastq_ch )
        adapter_ch          =   FIND_ADAPTERS   ( fastq_ch, fastqc_ch.fastqc_data )
        trimmed_fastq_ch    =   FASTP           ( fastq_ch, adapter_ch )
        
    emit:
        trimmed_fastq_ch.trimmed_fastq
}

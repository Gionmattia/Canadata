

process FASTQC {

    tag 'medium'

	publishDir "${params.output_dir}/fastqc", mode: 'copy'

    // errorStrategy  { task.attempt <= maxRetries  ? 'retry' :  'ignore' }
	
    // There's an issue...FASTQC can work with two files, but the output generation is messed up!
    // It tries to save both files together, which messes up everything.
    // "*_fastqc*/fastqc_data.txt", emit: fastqc_data
    // Instead it should process each singularly...
    // so...maybe change the tuple structure? I mean, we are processing two at the time...but where do I put the FASTQC report
    // so it is fed to the next step?


	input:
        tuple val(sample_id), file(fastq)

	output:
	    path "*_fastqc.html", emit: fastqc_html
        path "${fastq.baseName}_fastqc/fastqc_data.txt", emit: fastqc_data //should be this emitted as a tuple?

    script:
        """
        fastqc --extract -q $fastq --adapters ${projectDir}/scripts/adapter_list.tsv --dir ${projectDir}/data/temp
        """
}
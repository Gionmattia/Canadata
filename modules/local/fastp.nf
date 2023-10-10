process FASTP {
    publishDir "${params.output_dir}/fastp", mode: 'copy', pattern: '*.json'
    publishDir "${params.output_dir}/fastp", mode: 'copy', pattern: '*.html'

    // errorStrategy 'ignore'
    

    // Currently takes as input the tuple, with the paired files.
    input:
    tuple val(sample_id), file(raw_fastq)  // path(raw_fastq) 

    

    // Change the output so it gets the two files in output
    output:
    path ("${raw_fastq.baseName}_clipped.fastq"), emit: trimmed_fastq
    path '*1*.json', emit: fastp_json_1
    path '*1*.html', emit: fastp_html_1

	
	script:

    def reads1 = raw_fastq[0]
    def reads2 = raw_fastq[1]
    // this does not work.
    // NEXT APPROACH: I am just passing all the 4 arguments to the bash script
    // (the read1, the read2, the adapter_report1 and the adapter_report2)

	"""

    file1 = "${params.output_dir}/adapter_reports/${sample_id}__1.fastq.gz_adpater_report.fa"
    file2 = "${params.output_dir}/adapter_reports/${sample_id}__2.fastq.gz_adpater_report.fa"


    # Initialize variables
    a=""
    b=""

    # Loop through both files line by line
    awk '{
        if (\$0 !~ /^>/) {
            getline x < "'$file2'"
            if (x !~ /^>/) {
                # Assign to variables a and b
                a=\$0
                b=x

                print "Processing: " \$0 " " x

                fastp \
                -i $reads1 \
                -I $reads2 \
                -o ${reads1.baseName}_clipped.fastq \
                -O ${reads2.baseName}_clipped.fastq \
                --length_required 20 \
                --adapter_sequence $a  \
                --adapter_sequence_r2 $b  \
    -           --json ${reads1.baseName}.json \
                --html ${reads2.baseName}.html \
                --json ${reads1.baseName}.json \
                --html ${reads2.baseName}.html

            }
        } else {
            getline x < "'$file2'"
        }
    }' $file1


    """
}

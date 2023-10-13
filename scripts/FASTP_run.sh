#!/bin/bash
reads1=$1
reads2=$2
file1=$3
file2=$4

adap1=$(awk '!/^>/{print; exit}' $file1)
adap2=$(awk '!/^>/{print; exit}' $file2)

fastp \
    -i $reads1 \
    -I $reads2 \
    -o ${reads1%.*}_clipped.fastq \
    -O ${reads2%.*}_clipped.fastq \
    --length_required 20 \
    --adapter_sequence $adap1 \
    --adapter_sequence_r2 $adap2 \
    --json ${reads1%.*}.json \
    --html ${reads2%.*}.html

#!/bin/bash

# Check if both files are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: ./interleave_fasta.sh <file1.fasta> <file2.fasta>"
    exit 1
fi

# Variables to hold file names
file1=$1
file2=$2

# Check if files exist
if [ ! -f "$file1" ]; then
    echo "$file1 not found!"
    exit 1
fi

if [ ! -f "$file2" ]; then
    echo "$file2 not found!"
    exit 1
fi

# Initialize variables
a=""
b=""

# Loop through both files line by line
awk '{
    if ($0 !~ /^>/) {
        getline x < "'$file2'"
        if (x !~ /^>/) {
            # Assign to variables a and b
            a=$0
            b=x

            print "Processing: " $0 " " x

            # Run FASTP tool here with a and b
            # Replace this comment with your FASTP command, for example:
            # system("FASTP " a " " b)
        }
    } else {
        getline x < "'$file2'"
    }
}' $file1

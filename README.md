
# Project Canadata: the role of eIF4-BP in anoikis.

## Introduction 

The following pipeline allows the user to process polisome-profiling and RNA-seq data in order to obtain gene-count files.
It has been firstly used in the scope of "Project Canada".

[NB. Add the "DTA" part in the future...?]


## Requirements 
This pipeline is run using the following container methods
| Method      | Instructions                                                                                   |
| ----------- | ---------------------------------------------------------------------------------------------- |
| Singularity | [docs.syslabs.io](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)              |
| Conda       | [docs.conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  |

## Pipeline Structure

The main.nf pipeline comprises two different sub-workflows, named "preprocessing.nf" and "quantification.nf".
- Preprocessing.nf collects the raw files, performs quality checks and removes adapters and contaminants (more details in each process section).
- Quantification.nf takes PE files and proceeds to ontain transcript counts. These are then converted to gene counts.

## Required files for the analysis
This pipeline requires several files (annotations, indexes or similar) in order to be executed.
Given that many research lab could already have such files, the commands to obtain or create them have not been added to the main pipeline.
In case users do not have such files (or do not know what those are), all the commands used to retrieve/create them have been compiled into a set of instructions inside this documentation.

## Downstream DTA analysis
In the specific scope of "Project Canadata", the gene count files have also been merged and then used as input for  Differential Translation Analysis.
Considering the study-specific nature of the analysis (and the need for user input) these steps have not been added to the main pipeline.
Instead, the instruction to re-create the virtual environment (as well as all the executable scripts used) have been added to this directory to ensure reproducibility.

All the steps carried are also clearly documented, so that users can re-execute the analysis on their own machine.

<br>


## Setup
##### Input files needed

- Indexes and annotation files.
- parameters.yml file.
- nextflow.config file.
- experimental design file.
- Your data.

For details about each input, see sections below.

<br>

## Indexes and annotation files needed

Several processes in this pipeline will require to use indexes or annotation files (bowtie.nf, salmon.nf, tximport.nf). This section explains how to create each file and where to store them.
The following tutorial will use the annotations and indexes used in this study as an example.


### 1) Bowtie.nf

This module requires an index of the RNAs that the user wants to remove. Such index needs to be built from scratch.
You will need to use the following libraries and tools: wget, gunzip, python3 (libraries: argparse, Bio, SeqIo), bowtie.
**NB. If you do not have bowtie on your machine, you can create a conda environment with it (see end of this section)**

First, download the latest gencode fasta files containing the RNAs sequences.

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz
```

Then, unzip the file.

```
gunzip gencode.vM33.transcripts.fa.gz
```

Run the script "extract_rRNAs_mtRNAs.py" from the script folder of this Github repo. This script requires two arguments: the input file and the name for the output file.

```
python3 <path>/<to>/<your_local>/<canadata_repository>/scripts/extract_rRNAs_mtRNAs.py gencode.vM33.transcripts.fa gencode_rRNAs.fa
```

Once you have it, you need to index this file using bowtie. This command takes in input the .fa file generated the previous step and requires the user to specify the output name

```
bowtie-build  <gencode_rRNAs.fa> <ebwt_base>
```

Afterwards, the last thing is to move all the files (if they are not already) in a new directory. Within the "data" directory of this github repo, create a directory named "bowtie_indexes" (Nb. keep the name unchanged) and move all the files there.

```
mkdir ./bowtie_indexes
mv <location_of_ebtw_indexes>/*.ebwt <path>/<to>/<canadata>/data/bowtie_indexes
```

**OPTIONAL**

If you do not have bowtie on your machine, you can use the file bowtie.yml to create a temporary conda environment where to have it. Beware you will have to install the tools and libraries needed there, if you want to execute all the commands detailed above in that environment.

```
conda env create -f <path_to_canadata>\conda\bowtie.yml
```

_Congrats, you build the bowtie index!_


### 2) Salmon.nf

Salmon will require you to build an index for it to work. This code is basically the tutorial given by the salmon documentation, but broken down.

First, create a folder inside data to download the file into, then download the files needed: the genome.fa and the trasncripts.fa

**NB. If you executed the step above, you probably will aready have a .fa file for the transcripts. You can use that instead**

```
mkdir <path_to_canadata>/data/mouse_gencode_releases
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz
```

Then concatenate the two files into a "gentrome" file.

```
cat gencode.vM33.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > /data2/Canadata/Canadata/data/salmon_index/gentrome.fa.gz
```
Make a folder to store your salmon index:

```
mkdir <path_to_canadata>/data/salmon_index
```

Afterwards, you need to run salmon to build the index. You can do so by running salmon directly from the singularity container, following this code:

```
sudo singularity run --bind <path_to_canadata>/data/salmon_index <path_to_canadata>/singularity/salmon:1.10.1--h7e5ed60_0 salmon index -t <path_to_canadata>/data/salmon_index/gentrome.fa.gz -d <path_to_canadata>/data/salmon_index/decoys.txt -p 12 -i salmon_index --gencode
```

_Congrats, you build the salmon index!_

**NB2. The salmon index is built with a defo size of -k 31. This is also applied here. You might want to change that value depending on the length of your reads though.**


### 3) Tximport.nf

This module needs an annotation file for its .R script to work, which should be downloaded into the _data/mouse_gencode_releases_ subdirectory.

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.chr_patch_hapl_scaff.annotation.gtf.gz
```

_Congrats, you have made it!_



<br>


## Parameters.yml file

You can see an example of such file in the canadata repository. This file holds several variables (hard-coded paths) used by the nextlfow pipeline.
The pipeline also assumes you are running _main.nf_ from the canadata folder, once dloned locally.
- _input_dir_ -> The directory holding your fastq files.
- _output_dir_ -> The directory holding your output files. It is **mandatory** for the user to create the directory as a subdirectory of _data_ and to name it _output_
- _rRNA_index_ -> The path to the rRNA/mtRNA indexes you built in the previous section.
- _salmon_index_ -> The path to the salmon indexes you built in the previous setion.
- _singularity_path_ -> The path to the singularity image used to create the salmon container.
- _annotation_path_ -> The path to the annotation files (.gtf) that is needed by tximport.nf).

**NB. It is advisable to give the hard-coded path to the _salmon_index_ variable, since the machine could have issues when builded to bind the directory with the index.**


<br>

## Nextflow.config file

This file stores the configurations used by the pipeline when running. For reproducibility purposes, it should stay the same.
HOWEVER, this is also the file where the allocation of computational resources for each process are determined.
Before running the pipeline, check your capabilities of your local machine. You might want to adjust the following two parameters for each process:
- _maxForks_ -> Determines how many times the process can be run in parallel at the same time.
- _cpus_ (under the section profiles/standard/withLabel of the config file) -> It allocates by default 12 cpus to each process. You might want to reduce this number or come up with a different profile.

Each process is allowed _at most_ 12 CPUs per execution (some won't even get past 1, that is just a cap). Be aware that this value is _per process_: if you allow a process to have maxForks=2, then that process could end up suing up to 24 CPUs in total (12 per each of the two parallel executions).
Some processes will not go past 1CPU per execution (like FASTQC) but others will and could cause your system to crush!

In addition, once a file is ready to move to the next process, it will, thus increasing the resource consumptions!
For instance, let's assume you are processing 10 files. 6 could be in the first process (6 CPUs used), 2 could be in the second (12 CPus each) and the last 2 in the third (12 CPUs each again), thus determining a total consumption of 42 CPUs! Our server was fine, but that could be too much for other machines. If that's the case, adjust the value of _cpus_ from 12 to something more manageable by your system.

Be aware of your system capabilities before running this pipeline.


### How to create the bowtie index 
```
sudo singularity build singularity/pipeline Singularity
```
Then as the profile `singularity` specifies `container = 'singularity/pipeline'` use the following to execute:
```
nextflow run main.nf -profile singularity
```




##### Conda 
Create a conda definition yaml file [eg. here](conda/example.yml)
```
nextflow run main.nf -profile conda
```

## Usage
Call the pipeline directly
```
nextflow run main.nf
```

Run with all the frills
```
bash scripts/run-w-frills <params-file> <profile name from nextflow.config>
```
Example
```
bash scripts/run-w-frills example_parameters.yml standard
```


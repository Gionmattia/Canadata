
# Project Canadata: the role of eIF4-BP in anoikis.

## Introduction 

The following pipeline allows the user to process polisome-profiling and RNA-seq data in order to obtain gene-count files.
It has been firstly used in the scope of "Project Canada".

[NB. Add the "DTA" part in the future...?]

## Pipeline Structure

The main.nf pipeline comprises two different sub-workflows, named "preprocessing.nf" and "quantification.nf".
- Preprocessing.nf collects the raw files, performs quality checks and removes adapters and contaminants (more details in each process section).
- Quantification.nf takes PE files and proceeds to ontain transcript counts. These are then converted to gene counts.

## Downstream DTA analysis
In the specific scope of "Project Canadata", the gene count files have also been merged and then used as input for  Differential Translation Analysis.
Considering the study-specific nature of the analysis (and the need for user input) these steps have not been added to the main pipeline.
Instead, the instruction to re-create the virtual environment (as well as all the executable scripts used) have been added to this directory to ensure reproducibility.

All the steps carried are also clearly documented, so that users can re-execute the analysis on their own machine.

<br>

---
# Setup 


## Requirements 
This pipeline is run using the following container methods
| Method      | Instructions                                                                                   |
| ----------- | ---------------------------------------------------------------------------------------------- |
| Singularity | [docs.syslabs.io](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)              |
| Conda       | [docs.conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  |

## Input files needed
This pipeline requires several files (annotations, indexes or similar) in order to be executed.
Given that many research lab could already have such files, the commands to obtain or create them have not been added to the main pipeline.
In case users do not have such files (or do not know what those are), all the commands used to retrieve/create them have been compiled into a set of instructions inside this documentation.

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

**NB2. The salmon index is built with a defo size of -k 31. This is also applied here. You might want to change that value depending on the length of your reads though. This pipeline generates FASTQC reports which you can use to that end.**


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

<br> 


## Experimental Design file

This file is needed for the DTA analysis. It is a tab-separated .txt file that will allow the packages used in the DTA to correctly interpret your data.
**NB.** The "sampleID"s you need to insert at this point are the names of each sample you have in your analysis, as named at the end of the main.nf pipeline.

The first line fo the file must follow this structure (remember, the space between the name is a <tab>!):

```
SampleID	Condition	SeqType	Batch	Cell
```
Then, you just compile each row as if you were writing a tab-separated table. See below

```
SampleID	Condition	SeqType	Batch	Cell
Rep1_BPanoi_input_CACCTTAC-TAGTGCCA	anoi	RNA	1	BP
Rep1_BPanoi_poly_TTCTCGAC-ACCAAGCA	anoi	RIBO	1	BP
[...]
```
Remember this example is specific to the context of the Canadata Project, were "SeqType", "Cell", and "Condition" were contrasts in the analysis. 
If you need more information on this file for your own experiments, look at the DeltaTE documentation at (https://github.com/SGDDNB/translational_regulation)

<br>

---

# Pipeline Usage
Call the pipeline directly from the cloned repository. Remember to specify the usage of the params_files.yml

```
nextflow run main.nf -params-file example_parameters.yml
```

The pipeline will run across all the steps. Final output will be quantification files (gene_counts.csv files).
In order to execute the DTA analysis, it will be necessary to run these files through some extra scripts.
Given the study-specific nature of these analyses, such passages have not been added to the pipeline, but all the steps and scripts are still included in this github repo for reproducibility (see section below).

<br>

---

# DTA analysis

This section details all the scripts used for the DTA analysis, as well as some custom scripts used to generate graphical or tabular outputs for consulation.

## Merging the gene count files

First thing first, you will need to combine the samples together by SeqType: you want a file with all your gene counts from RNA-seq data and another for the ones from Rivo-seq/polysome profiling.

This can be done with the script "tsv_merge.py" in the _scripts_ folder of this repository. To run the script you will need python3 and the following libraries:
- pandas
- argparse
- os

You can then execute the script from the terminal as it follows:

```
python3 ./scripts/tsv_merge.py <input_directory> <keyword> <output>
```
- _input_directory_ is basically the directory where you have all your gene counts after the TXIMPORT module has converted them. It should be the "./data/output/gene_counts" if you execute the pipeline as it is.
- _keyword_ is the substring within the name of each gene count file that determines whether it is from RNAseq or RIBO-seq/Poylsome profiling (in this study, it was "input" and "poly", respectively)
- _output_ is just the name of the output file. You can also specify a directory for it.

At the end of this step, you should have two files in .tsv format. One containing all the gene counts (rows) for each sample (columns) for RNA-seq samples, the other for Ribo-Seq/Polysome profiling.

## Running the DTA Analysis

We are finally at the end, the DTA analysis. You will need to set up a conda environment to hold all the packages you need. There's a .yml file for it in the _conda_ subdirectory (_delthena.yml_).
You can set everything by using this command in the terminal:

```
conda env create -f ./conda/delthena.yml 
```

Now just activate the environment and execute the scripts for the analysis. 

```
conda activate delthena
```

Again, these scripts are specific to this study. If you want to use them for your own analysis, you will need to adapt them to your study (specifically, your contrasts and your choice of baseline for each).

```
Rscript ./scripts/DeltaTE_ORFik_Analysis(by_cell).R <path_to_RNA-seq_count_file> <path_to_RIBO/Polysome_count_file> <path_to_sample_info.txt> <Batch_correction_T_or_F>
Rscript ./scripts/DeltaTE_ORFik_Analsysis(by_condition).R <path_to_RNA-seq_count_file> <path_to_RIBO/Polysome_count_file> <path_to_sample_info.txt> <Batch_correction_T_or_F>
```

These scripts need 4 arguments, which are in order:
- The path to the RNA-seq gene count file
- The path to the Ribo-seq/Polysome profiling gene count file
- The path to your *experimental design file* (see section above), here called _sample_info.txt_
- A boolean (T/F) which will inform the script whether a batch analysis is necessary or not. You can see if a batch effect is present in your data by performing a PCA analysis (more info below).

_Congrats, you completed the DTA analysis and the pipeline so far!_

<br>

---

# DLC Content

These are things not really required for a DTA analysis pipeline, but still added for future reference 

Scripts to perform a PCA analsys
Scripts to plot the data and convert it to gene names/symbols
[add later]





Run with all the frills
```
bash scripts/run-w-frills <params-file> <profile name from nextflow.config>
```
Example
```
bash scripts/run-w-frills example_parameters.yml standard
```


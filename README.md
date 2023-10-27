
# Project Canadata: the role of eIF4-BP in anoikis.

## Introduction 

The following pipeline allows the user to process polisome-profiling and RNA-seq data in order to obtain gene-count files.
It has been firstly used in the scope of "Project Canada".

[NB. Add the "DTA" part in the future...?]


## Requirements 
This pipeline can be run using each of the following container methods
| Method      | Instructions                                                                                   |
| ----------- | ---------------------------------------------------------------------------------------------- |
| Singularity | [docs.syslabs.io](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)              |
| Conda       | [docs.conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  |

## Pipeline Structure

The main.nf pipeline comprises two different sub-workflows, named "preprocessing.nf" and "quantification.nf".
- Preprocessing.nf collects the raw files, performs quality checks and removes adapters and contaminants (more details in each process section).
- Quantification.nf takes PE files and proceeds to ontain transcript counts. These are then converted to gene counts.

### Extra appendages

#### Required files for the analysis
This pipeline requires several files (annotations, indexes or similar) in order to be executed.
Given that many research lab could already have such files, the commands to obtain or create them have not been added to the main pipeline.
In case users do not have such files (or do not know what those are), all the commands used to retrieve/create them have been compiled into a set of instructions inside this documentation.

#### Downstream DTA analysis
In the specific scope of "Project Canadata", the gene count files have also been merged and then used as input for  Differential Translation Analysis.
Considering the study-specific nature of the analysis (and the need for user input) these steps have not been added to the main pipeline.
Instead, the instruction to re-create the virtual environment (as well as all the executable scripts used) have been added to this directory to ensure reproducibility.

All the steps carried are also clearly documented, so that users can re-execute the analysis on their own machine.


##### Input files needed


## Setup
##### Input files needed

- Indexes and how to create them (step by step guide).
- parameters.yml file.
- nextflow.config file.
- -experimental design file.
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


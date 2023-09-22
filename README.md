# Nextflow-Template
Template repository for setting up new Nextflow pipelines

This template has been informed by many other Nextflow pipelines and templates. Some are much more advanced than this template. Please see:
- [Awesome-Nextflow](https://github.com/nextflow-io/awesome-nextflow)
- [Australian BioCommons Template](https://github.com/AustralianBioCommons/Nextflow_DSL2_template)
- [epi2me labs](https://github.com/epi2me-labs/wf-template)
- [FredHutch](https://github.com/FredHutch/workflow-template-nextflow)
- [And Many More](https://github.com/search?q=Nextflow%20Template&type=repositories)

## General Setup Notes

#### Tower 
tower.nf allows GUI tracking of workflow progress
Set it up using the instructions [here](https://help.tower.nf/22.4/getting-started/usage/#nextflow-with-tower)

<br>
<br>
<br>

# ||=====================================||
# ||   Template README               
# ||   <b>Delete Everything Above</b>
# ||=====================================||

# <PIPELINE NAME>
## Introduction 

[Describe here what this pipeline does]

## Requirements 
This pipeline can be run using each of the following container methods
| Method      | Instructions                                                                                   |
| ----------- | ---------------------------------------------------------------------------------------------- |
| Singularity | [docs.syslabs.io](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)              |
| Docker      | [docs.docker.com](https://docs.docker.com/engine/install/)                                     |
| Conda       | [docs.conda.io](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)  |


## Setup
##### Singularity
```
sudo singularity build singularity/pipeline Singularity
```
Then as the profile `singularity` specifies `container = 'singularity/pipeline'` use the following to execute:
```
nextflow run main.nf -profile singularity
```

##### Docker
```
docker build . -t pipeline-image
```
Then as the profile `docker` specifies `container = 'pipeline-image:latest'` use the following to execute:
```
nextflow run main.nf -profile docker
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


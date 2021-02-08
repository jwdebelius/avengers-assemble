# Avengers Assemble

This repository contains the analysis workflow for Debelius et al (doi: []()). 


This 

## Running the pipeline

### Installation

The analysis was done based on a qiime2-2020.11 enviroment with [RESRIPt]() and [q2-sidle]() installed. To build this enviroment, install [qiime2-2020.11]() according to your system requirements. Then, follow the [q2-sidle]() installation instructions. You will also need snakemake > 5.3.

```bash
conda install --channel bioconda snakemake=5.3
```


### Input Files


...

### Running the workflow

#### Workflow steps

##### Simulation

1. **Prepare the Silva 128 database using RESCRIPt.** This set of commands downloads the database, filters degenerate sequences, Eukeryotes, and sequences with undefined phyla. 


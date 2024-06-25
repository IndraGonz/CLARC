# Connected Linkage and Alignment Redefinition of COGs (CLARC)

A tool that uses sequence identity, within-population linkage patterns and functional annotations to identify and reduce the over-splitting of genes into multiple clusters of orthologous genes (COGs) in a pangenome analysis. 

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quickstart](#quickstart)
   - [Basic usage](#basic-usage)
   - [Test data](#test-data)
4. [Pipeline workflow description](#pipeline-workflow-description)
5. [Full usage](#full-usage)
6. [Outputs](#outputs)
7. [Raising issues](#raising-issues)
8. [Citing](#citing)

## Introduction

🚧🚧 Introduction coming soon 🚧🚧

## Installation

### Requirements

For now, CLARC can only be installed from source by copying this github repository. Also, all dependencies are installed using a conda environment. So, the installation assumes that you have git and conda installed. 

For quick and simple instructions of how to locally install conda with miniconda follow [this link](https://docs.anaconda.com/free/miniconda/#quick-command-line-install). 

Then you can simply install git with:

```bash
conda install git
```

### Install from source

Clone repository from GitHub:

```bash
git clone https://github.com/IndraGonz/CLARC.git
```

The folder will be copied into your working directory, so make sure to navigate to your desired path for the installation.

Now, navigate to the cloned folder:

```bash
cd CLARC
```

Create conda environment with the necessary dependencies:

```bash
conda env create --file envs/clarc_env.yml
```

Activate environment:

```bash
conda activate clarc_env
```

Now from the active clarc environment, in the CLARC folder, run the installation command:

```bash
python setup.py install
```

That should do it, now you should have your very own CLARC installation within the ```clarc_env``` conda environment!

To confirm that the installation worked type:

```bash
clarc --version
```

If the version appears, then the installation was successful!

## Quickstart

### Basic usage

#### Description of inputs

CLARC's default is to take the results of a [Roary](https://github.com/sanger-pathogens/Roary) pangenome analysis as input. However, it can also accept results from the [Panaroo](https://github.com/gtonkinhill/panaroo) pipeline if the ```--panaroo``` flag is specified. More information on the panaroo usage in the [full usage](#full-usage) section of this documentation. Basic usage will assume that Roary was used to generate the pangenome.

From Roary, CLARC uses two output files: The ```gene_presence_absence.csv``` which contains a presence absence matrix with all the COGs identified by Roary and ```pan_genome_reference.fa``` which is a fasta file containing the representative sequences of these genes. The representative sequence is the longest instance of that COG across all samples where that COG was called. The .csv file contains ALL COGs called by Roary and not only core/accessory. Do not rename these files, as CLARC will look for them by name.

Additionally, CLARC needs a text file with the names of the samples in your population of interest. This file should be named ```needed_sample_names.txt```. The 'accessory gene' definitions and the linkage constraints will be calculated based on this population of samples. 

For example, if you ran the pangenome analysis using samples from various geographic locations/datasets, you can feed CLARC only the samples in a particular geographic location you are interested in.

Nonetheless, if you did not use samples from multiple distinct populations in your original pangenome analysis, you can just include the names of all the samples used to build the pangenome in the ```needed_sample_names.txt``` input file. 

Make sure these 3 input files are in the same folder, since you will specify the input folder path when running CLARC.

#### Command

Now the clarc command can be run from the terminal (within the clarc_env environment):

```bash
clarc --input_dir {path to folder with input data} --output_dir {path to folder where the clarc results will be stored}
```

### Test data

Instructions of how to run CLARC on a given set of input files can be found in the [pangenome analysis](https://github.com/IndraGonz/nfds-tutorial?tab=readme-ov-file#pangenome-analysis) section of this broader tutorial. Specifically, the CLARC exercise is in the *Running CLARC* subsection. 

In this example, CLARC uses Roary results from a multi-population pangenome built with >8000 S. pneumoniae samples collected across 7 different locations. It specifically refines the COG definitions for one of the datasets, which was collected in the Southwest of the United States.

## Pipeline workflow description

🚧🚧 All descriptions coming soon 🚧🚧

### Creating necessary inputs

#### Accessory and core gene filtering

#### Linkage matrix calculation

#### All vs. all BLAST for accessory COGs

#### Functional annotation with the EggNOG database

### CLARC algorithm

## Full usage

```bash
Usage:   clarc --input_dir {input data path} --output_dir {output path} [options]

Options: --acc_upper INT   Upper frequency threshold for accessory gene filtering, if not specified default is <= 0.95 (95%)
         --acc_lower INT   Lower frequency threshold for accessory gene filtering, if not specified default is > 0.05 (5%)
         --core_lower INT   Lower frequency threshold for core gene filtering, if not specified default is > 0.95 (95%)
         --panaroo   Flag specifying that panaroo inputs will be used, if not provided the program will assume the pangenome inputs are from Roary
         --filter-only   Flag to only run the initial accessory/core COGs filtering step, not the CLARC clustering algorithm. 
         --options   Display these options
         --version   Display CLARC version currently installed
```

### Inputs for Panaroo

## Outputs

### ```clarc_results``` folder (main outputs)

#### ```clarc_cluster_summary.txt```

#### ```accessory_cluster_cogs.txt``` and ```core_cluster_cogs.txt```

#### ```accessory_cluster_summary.csv``` and ```core_cluster_summary.csv```

#### ```clarc_condensed_presence_absence.csv```

#### ```clarc_acc_cog_seqs.fasta```

Additional output files are:

### ```population_accessory_presence_absence.csv```

csv file containing the presence absence matrix for the _accessory_ genes within the samples specified in the ```needed_sample_names.txt``` input file. So, only COGs present in between 5-95% of those samples will be included here (if run on default parameters). This is just a filtered output from the original pangenome analysis, before running the CLARC clustering algorithm.

### ```accessory_rep_seqs.fasta```

fasta file containing the representative DNA sequences of all _accessory_ COGs identified in the given population of interest. The entries in this file will correspond to the COGs in the ```population_accessory_presence_absence.csv``` output file. The COG name (as given by the pangenome tool) will be the fasta entry identifier. 

### ```accessory_rep_protein_seqs.fasta```

fasta file containing the representative protein sequences of all accessory COGs identified in the given population of interest. It is just the entried in the ```accessory_rep_seqs.fasta``` file translated from DNA to protein sequences. The pipeline will internally use these protein sequences to functionally annotate the accessory COGs using the EggNOG database.

### ```population_core_presence_absence.csv```

csv file containing the presence absence matrix for the _core_ genes within the samples specified in the ```needed_sample_names.txt``` input file. So, only COGs present in over 95% of those samples will be included here (if run on default parameters). This is just a filtered output from the original pangenome analysis, before running the CLARC clustering algorithm. Essentially it it the ```population_accessory_presence_absence.csv``` output, but for core genes instead of accessory.

### core_rep_seqs.fasta

Like, the ```accessory_rep_seqs.fasta``` output file, this is a fasta file containing the representative DNA sequences of all _core_ COGs identified in the given population of interest.



### population_accessory_presence_absence.csv

### accessory_rep_seqs.fasta

### accessory_rep_protein_seqs.fasta

### population_core_presence_absence.csv

### core_rep_seqs.fasta

### linkage folder

### acc_blastn folder

### eggnog folder




## Raising issues

If you have trouble with any of the steps shown here, please let me know! You can email me at igonzalezojeda@g.harvard.edu or submit an issue through GitHub. Also feel free to reach out with any general questions/suggestions!

## Citing

🚧🚧 Coming soon 🚧🚧


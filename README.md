# CLARC - Connected Linkage and Alignment Redefinition of COGs

<!-- Using HTML for more control over the size -->
<img src="https://github.com/IndraGonz/CLARC/blob/main/CLARC.egg-info/CLARC_logo.svg" alt="CLARC logo" width="450" height="270">

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

```
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

Running clarc will create various subdirectories containing the intermediate files necessary to run the clustering algorithm. These include the linkage matrices created, the eggnog functional annotations of the original accessory COGs and the results of the all vs all BLAST of these accessory COGs. 

However, the final results are located in the subdirectory named ```clarc_results```. In the future I will add an option to suppress the intermediate outputs, and I will update this document when I do. 

With that said, let's start by reviewing the final output files of the pipeline.

### clarc_results folder (final output files)

The ```clarc_results``` folder will contain the following files:

   - clarc_cluster_summary.txt

      simple txt file containing a summary of the 'same gene' clusters found by CLARC, if any. It is in this format:

      ```
      Core COG clusters: X
      Unique COGs in core clusters: X
      Accessory COG clusters: X
      Unique COGs in accessory clusters: X
      ```

      If CLARC finds no clusters, then this file will just contain a message indicating that no clusters were found.
      
   - clarc_condensed_presence_absence.csv

     csv file containing the presence absence matrix for the CLARC re-defined accessory COG definitions. COGs in core clusters are removed and those in accessory clusters are collapsed into a        new cluster. The name of the new clusters will be the names of the original COGs forming the cluster, joined by slashes. For example: 'COG1/COG2/COG3/' for a 3 COG cluster identified by         CLARC.
     
   - clarc_acc_cog_seqs.fasta

      fasta file containing the representative sequences of the new CLARC accessory COG definitions. The entries in this file will correspond to the COGs in the                                        ```clarc_condensed_presence_absence.csv``` output file. The COG name will be the fasta entry identifier, with new clusters having the appropiate cluster name. For clusters, the longest          sequence between the individual representative sequences of each COG is chosen as the cluster representative sequence.
     
   - accessory_cluster_cogs.txt and core_cluster_cogs.txt

     txt files that include a list of the unique original COGs present in the accessory and core clusters identified by CLARC. No information on the actual clusters (which COGs are connected to      which) is included here, just the list of the original COGs that belong to any cluster.

   - accessory_cluster_summary.csv and core_cluster_summary.csv

   csv file specifying the compositions of each accessory and core cluster identified by CLARC. The number of columns will depend on the number of members in the largest cluster. The original      frequency in the population of each COG member is also included (freq_cogX) as well as the new frequency of the re-defined cluster (freq_sum).
   
Intermediate output files are:

In the filtering step (before CLARC re-definition), 5 files are created:

   - population_accessory_presence_absence.csv

      csv file containing the presence absence matrix for the _accessory_ genes within the samples specified in the ```needed_sample_names.txt``` input file. So, only COGs present in between          5-95% of those samples will be included here (if run on the default parameters). This is just a filtered output from the original pangenome analysis, before running the CLARC clustering         algorithm.
     
   - accessory_rep_seqs.fasta

      fasta file containing the representative DNA sequences of all _accessory_ COGs identified in the given population of interest. The entries in this file will correspond to the COGs in the        ```population_accessory_presence_absence.csv``` output file. The COG name (as given by the pangenome tool) will be the fasta entry identifier. 

   - accessory_rep_protein_seqs.fasta

      fasta file containing the representative protein sequences of all accessory COGs identified in the given population of interest. It is just the sequences in the                                  ```accessory_rep_seqs.fasta``` file translated from DNA to protein. The pipeline will internally use these protein sequences to functionally annotate the accessory COGs using the EggNOG         database.

   - population_core_presence_absence.csv

      csv file containing the presence absence matrix for the _core_ genes within the samples specified in the ```needed_sample_names.txt``` input file. So, only COGs present in over 95% of           those samples will be included here (if run on default parameters). This is just a filtered output from the original pangenome analysis, before running the CLARC clustering algorithm.           Essentially it is the ```population_accessory_presence_absence.csv``` output, but for core genes instead of accessory.

   - core_rep_seqs.fasta

      fasta file containing the representative DNA sequences of all _core_ COGs identified in the given population of interest. Likewise, this file is equivalent to the                                ```accessory_rep_seqs.fasta``` output file, but for core genes instead of accessory.

### linkage folder

   - acc_pXX_matrix.csv

     csv files that are pairwise matrices with counts of genomes that show the different states of co-occurrence between every possible pair of accessory COGs. The possible states are: genomes       where both COGs are observed together (P11), genomes where neither COG is observed (P00), or genomes where only one of the COGs is observed (P10 and P01). Thus, 4 of these matrices are          generated, one per co-occurrence state.

   - acc_linkage_co-occur.csv

     csv file containing a pairwise metric of linkage. It is internally calculated using the counts in the previous ```acc_pXX_matrix.csv``` matrices, with the following equation:

     <img width="206" alt="image" src="https://github.com/IndraGonz/CLARC/assets/70214497/49f8b26f-f9b9-4463-b12f-92c7d1895f58">

     This provides a normalized linkage metric where:
     
        - L > 0 implies a positive correlation between the genes in the pair (co-occur often)
        - L < 0 implies a negative correlation between the genes in the pair (exclude each other)
        - L = 0 implies no correlation between the genes in the pair

### acc_blastn folder

Various files are created here, since a BLAST database is built from the accessory cog sequence fasta file. However the only relevant file is the one containing the results for the all vs. all accessory COG BLAST:

   - blastn_acccogs_allvall.tsv

     tsv file of the BLAST results with 12 columns, with the following headers:

      ```python
     ['query_seq_ID', 'subject_seq_ID', 'percentage_indentical_matches', 'align_length','num_mismatches','gap_open','align_start_query','align_end_query','align_start_subject', 'align_end_subject','e-value','bit_score']
      ```
      
### eggnog folder

   - acc_cog_eggnog_annotations.csv

     csv file containing three columns that include the inputs and results of the EggNOG functional annotation. The first column in the accessory COG name, the second is the protein sequence         that was queried, and the third column is the EggNOG functional group identified for that accessory COG.

     Only the abbreviation for each COG functional group is shown, so here is a table including the meaning of those abbreviations:

     | Abbreviation | Function |
     |--------------|----------|
     | A | RNA processing and modification |
     | B | Chromatin Structure and dynamics |
     | C | Energy production and conversion |
     | D | Cell cycle control and mitosis |
     | E | Amino Acid metabolis and transport |
     | F | Nucleotide metabolism and transport |
     | G | Carbohydrate metabolism and transport |
     | H | Coenzyme metabolism |
     | I | Lipid metabolism |
     | J | Translation |
     | K | Transcription |
     | L | Replication and repair |
     | M | Cell wall/membrane/envelop biogenesis |
     | N | Cell motility |
     | O | Post-translational modification, protein turnover, chaperone functions |
     | P | Inorganic ion transport and metabolism |
     | Q | Secondary Structure |
     | T | Signal Transduction |
     | U | Intracellular trafficing and secretion |
     | V | Defense mechanisms |
     | Y | Nuclear structure |
     | Z | Cytoskeleton |
     | R | General functional prediction only |
     | S | Function unknown |
     | NA | No annotation - Protein had a hit in the EggNOG database, but the function is not annotated |
     | NH | No hit - Protein had no hit in the EggNOG database |
     | Mixed | Hit with multiple functions found |

     More information about these functional group classifications can be found in the [COG Database NIH website](https://www.ncbi.nlm.nih.gov/research/cog/#).

   - eggnog_group_cog_names.csv

     csv containing a summary of the COGs found in each category. The column headers are the categories, and each column has the list of COGs identified in that category. 

These functional annotation files (and all other intermediate files) are created for the accessory genes filtered from the original pangenome analysis results. However, one of the conditions that CLARC uses to create the same gene clusters is that the members of the cluster must have the same functional group. So, using these annotations the user can also obtain the functional annotated for the new CLARC accessory gene definitions.

## Raising issues

If you have trouble with any of the steps shown here, please let me know! You can email me at igonzalezojeda@g.harvard.edu or submit an issue through GitHub. Also feel free to reach out with any general questions/suggestions!

## Citing

🚧🚧 Coming soon, send me good vibes if you can spare them 🚧🚧


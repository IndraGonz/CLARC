# CLARC - Connected Linkage and Alignment Redefinition of COGs

<!-- Using HTML for more control over the size -->
<img src="https://github.com/IndraGonz/CLARC/blob/main/CLARC.egg-info/CLARC_logo.svg" alt="CLARC logo" width="450" height="270">

[![Python](https://img.shields.io/badge/Python-3.6-blue.svg)](https://www.python.org/)
[![MIT License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A tool that uses sequence identity, linkage patterns and functional annotations to identify and reduce the over-splitting of accessory genes into multiple clusters of orthologous genes (COGs) in a pangenome analysis. 

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

Current tools that infer bacterial pangenomes generally cluster annotated coding sequences into orthologous groups (COGs) to generate a gene presence absence matrix for the given population of genomes. However, strict identity cutoffs often inflate accessory gene estimates by misclassifying variants of a gene into separate orthologous groups. These misclassifications can significantly impact downstream analyses, especially those that rely on gene frequency estimates.  CLARC helps correct this over splitting of accessory genes into multiple orthologous groups, by identifying and condensing redundant COGs.

In its first step, CLARC identifies ‘same gene’ pairs by looking for genes that never co-occur in the same isolate and also share sequence identity and functional annotation. After this step, CLARC builds a graph where fully connected clusters represent gene variants that were erroneously split into different COGs. Finally, genes in these clusters are condensed to generate a refined gene presence absence matrix for the population. For a more detailed breakdown of the technical workflow see the [Pipeline workflow description](#pipeline-workflow-description) section in this repository.

In summary, CLARC is meant to compliment existing bacterial pangenome tools by polishing their COG definitions. As input, the pipeline currently takes the presence absence matrix generated with [Roary](https://sanger-pathogens.github.io/Roary/) (but can also accept input from [Panaroo](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02090-4)). We believe CLARC is particularly helpful for researchers that plan to perform downstream analyses that rely on COG frequencies, such as studying the evolutionary dynamics of accessory genes or running a panGWAS. 

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

CLARC's default is to take the results of a [Roary](https://github.com/sanger-pathogens/Roary) pangenome analysis as input. However, it can also accept results from the [Panaroo](https://github.com/gtonkinhill/panaroo) pipeline if the ```--panaroo``` flag is specified. More information on the panaroo usage in the [Inputs for Panaroo](#inputs-for-panaroo) subsection of this documentation.

From Roary, CLARC uses two output files: The ```gene_presence_absence.csv``` file which contains a presence absence matrix with all the COGs identified by Roary and ```pan_genome_reference.fa``` which is a fasta file containing the representative sequences of these genes. The representative sequence is the longest instance of that COG across all samples where that COG was called. The .csv file contains ALL COGs called by Roary and not only core/accessory. Do not rename these files, as CLARC will look for them by name.

Additionally, CLARC needs a text file with the names of the samples in your population of interest. This file should be named ```needed_sample_names.txt```. The 'accessory gene' definitions and the linkage constraints will be calculated based on this population of samples. 

For example, if you ran the pangenome analysis using samples from various geographic locations/datasets, you can feed CLARC only the samples in a particular geographic location you are interested in.

Nonetheless, if you did not use samples from multiple distinct populations in your original pangenome analysis, you can just include the names of all the samples used to build the pangenome in the ```needed_sample_names.txt``` input file. 

Make sure these 3 input files are in the same folder, since you will specify the input folder path when running CLARC.

#### Command

Now the clarc command can be run from the terminal (within the clarc_env environment):

```bash
clarc --input_dir {path to folder with input data} --output_dir {path to folder where the clarc results will be stored}
```

### Inputs for Panaroo

Because of the way in which Panaroo saves the representative sequences of each COG, the inputs will be slightly different. CLARC will need 3 files that are outputs from Panaroo: (1) The ```gene_presence_absence_roary.csv``` file that contains the presence absence matrix with all the COGs identified by Panaroo, in the Roary output format (2) The ```pan_genome_reference.fa``` file containing most COG representative sequences and (3) The ```gene_data.csv``` file that contains more detailed sequence information for all COGs called by Panaroo. This third output is necessary because Panaroo does not include all COG sequences in the ```pan_genome_reference.fa``` , and so an internal CLARC script uses the detailed information to generate a fasta file with all sequences.

As previously described with Roary, CLARC will need a text file with the names of the samples in your population of interest. This file should be named ```needed_sample_names.txt```. The 'accessory gene' definitions and the linkage constraints will be calculated based on this population of samples. 

Make sure these 4 input files are in the same folder, since you will specify the input folder path when running CLARC.

Now the clarc command can be run from the terminal, specifying the panaroo flag (within the clarc_env environment):

```bash
clarc --input_dir {path to folder with input data} --output_dir {path to folder where the clarc results will be stored} --panaroo
```

### Test data

Instructions of how to run CLARC on a given set of input files can be found in the [pangenome analysis](https://github.com/IndraGonz/nfds-tutorial?tab=readme-ov-file#pangenome-analysis) section of this broader tutorial. Specifically, the CLARC exercise is in the *Running CLARC* subsection. 

In this example, CLARC uses Roary results from a multi-population pangenome built with >8000 S. pneumoniae samples collected across 7 different locations. It specifically refines the COG definitions for one of the datasets, which was collected in the Southwest of the United States.

## Pipeline workflow description

CLARC uses a custom clustering algorithm to identify and reduce redundancy in COG definitions. The first step in the clustering process is to identify accessory COG pairs that appear to be the ‘same unit’ using 3 constraints: (1) COGs never co-occur in the same isolate, (2) COGs meet a custom cutoff of nucleotide sequence similarity and, (3) COGs get classified in the same EggNOG functional group. 

After this step, the algorithm builds a graph where each node represents a COG, and two nodes are connected by an edge if the constraints were met for that COG pair in the previous step. All three of the constraints must be met for two COGs to be connected by an edge in this graph. The algorithm then looks for fully connected clusters of 2 or more COGs that appear to be the same unit. By fully connected, we mean clusters where all COGs have connections to all other COGs in the cluster. This ensures that all fully connected clusters are mutually exclusive, and helps prevent false positives for low frequency COGs. COGs that are not in fully connected clusters within the graph are not modified further. Finally, the pipeline redefines the COGs by summing the presence absence matrix for all COGs within previously identified ‘same unit’ clusters. 

<img src="https://github.com/IndraGonz/CLARC/blob/main/CLARC.egg-info/github_workflow_1.png" alt="CLARC algorithm workflow" width="900">

This clustering algorithm is implemented in the CLARC bioinformatic tool. When using the tool, the user inputs a csv file with the COG-presence-absence matrix previously generated by a pangenome tool, along with the fasta file containing the COG representative sequences for that pangenome analysis (which is also an output from current pangenome tools). The current version of CLARC (v.1.1.2) supports raw inputs from Roary and Panaroo. However, the results from any pangenome analysis can be run with CLARC, as long as the inputs are formatted like the results from Roary or Panaroo. In addition to these two files, CLARC requires the user to provide a text file with a sample ID list for the genomes that compose the population of interest. This specified population can simply be all samples used to build the pangenome analysis, or a subset of samples if a distinct population within the whole set is known (e.g. samples collected in a specific geographic location).

<img src="https://github.com/IndraGonz/CLARC/blob/main/CLARC.egg-info/github_workflow_2.png" alt="CLARC tool workflow" width="900">

### Creating necessary inputs

#### Accessory and core gene filtering

An internal python script (filtering_acc_core.py/filtering_acc_core_panaroo.py) filters the original pangenome presence absence matrix to generate filtered presence absence csv files of the accessory and core genomes. This is done by calculating the frequency of each COG in the specified population of samples and then filtering for core and accessory genes using user-specified frequency thresholds. If run on the default parameters, CLARC will define accessory genes as those in between 5-95% frequency and core genes as those with >95% frequency. The csv file with the accessory gene presence absence matrix in this step will be used to calculate the internal CLARC algorithm inputs. 

#### Linkage matrix calculation

In this step, a python script (get_linkage_matrix.py) counts the different states of co-occurrence between every possible pair of accessory COGs. The possible states are: genomes where both COGs are observed together (P11), genomes where neither COG is observed (P00), or genomes where only one of the COGs is observed (P10 and P01). Thus, 4 of these matrices are generated, one per co-occurrence state.

#### All vs. all BLAST for accessory COGs

A bash script (acccog_blastn.sh) uses all accessory COG representative sequences to build a BLASTN database and perform an all vs. all BLAST comparison. BLASTN version 2.15.0 is part of the internal conda environment created with the initial CLARC build during installation. 

#### Functional annotation with the EggNOG database

CLARC uses a custom python script (eggnog_annotations.py) that serves as a wrapper of the [EggNOG v.5.0 database](http://eggnog5.embl.de/#/app/home). Accessory COG representative sequences are translated into their protein sequences and submitted as queries to the database [query interface](http://eggnog5.embl.de/#/app/seqscan). This allows for the remote functional annotation of COGs without having to download the EggNOG database locally (as is done when using tools like [EggNOG Mapper](http://eggnog-mapper.embl.de/)).

### CLARC algorithm

Finally, a python script (clarc_condense.py) containing the CLARC clustering algorithm uses the internal inputs generated to condense the COG definitions and create the outputs for the re-defined COG list. 

## Full usage

```
Usage:   clarc --input_dir {input data path} --output_dir {output path} [options]

Options: -h, --help  Show help message with options
         --acc_upper   Upper frequency threshold for accessory gene filtering, if not specified default is <= 0.95 (95%)
         --acc_lower   Lower frequency threshold for accessory gene filtering, if not specified default is > 0.05 (5%)
         --core_lower  Lower frequency threshold for core gene filtering, if not specified default is > 0.95 (95%)
         --ci, --clarc-identity BLASTN identity threshold CLARC uses as constraint to identity same gene clusters. Number from 0-100, default is 95%
         --panaroo   Flag specifying that panaroo inputs will be used, if not provided the program will assume the pangenome inputs are from Roary
         --filter-only   Flag to only run the initial accessory/core COGs filtering step, not the CLARC clustering algorithm. 
         --options   Display these options
         --version   Display CLARC version currently installed
         --dif, --delete-intermediate-files Delete intermediate files
         --merge MERGE [MERGE ...] Paths to CLARC results that will be merged, in the case of runs with multiple defined populations (folders where clarc was run that include the original data and clarc output subfolders)
```

## Outputs

Running clarc will create various subdirectories containing the intermediate files necessary to run the clustering algorithm. These include the linkage matrices created, the eggnog functional annotations of the original accessory COGs and the results of the all vs all BLAST of these accessory COGs. 

However, the final results are located in the subdirectory named ```clarc_results```. All the intermediate files can be cumbersome, and so you can suppress the output to only get the ```clarc_results``` folder by using the ```--dif``` or ```--delete-intermediate-files``` flag in the clarc command. 

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

   - ```roary_clarc_gene_presence_absence.csv```

     csv file that contains the modified version of the input pangenome csv presence absence file. This is in the format of the Roary csv output, since Panaroo also uses this format. The new         names of condensed clusters will be the old names joined by dashes. For example: 'COG1-COG2-COG3-' for a 3 COG cluster identified by CLARC.

     For each cluster, the other metadata information within the table (e.g. Annotation, Genome fragment, etc.) will be the entries of the original COGs, in order, joined by a comma. For             example, the 'Genome Fragment' column of the previous example 3 COG cluster would look like 'num1, num2, num3'. 

   - ```clarc_condensed_gene_pres_abs_binary.csv```

     csv file containing the 'binary version' of the condensed roary_clarc_gene_presence_absence.csv results. This means that instead of being formatted as the Roary/Panaroo output, it will be       formatted as a simple presence absence table where each row is a sample used in the pangenome analysis and each column is a COG. If the COG is present in the sample, the value                   will be a 1 and otherwise it will be 0.

   - ```clarc_population_acc_presence_absence.csv``` and ```clarc_population_core_presence_absence.csv```

      In addition to the full condensed matrix, CLARC will filter for the accessory and core genes only within the specified population (the samples in the ```needed_sample_names.txt``` input         file) and provide a csv file with the binary presence absence table for only these accessory/core genes. The frequency thresholds use to filter will be 5-95% for accessory genes and >95%        for core genes, unless different thresholds are specified when running clarc.
        
   - ```clarc_acc_cog_seqs.fasta``` and ```clarc_core_cog_seqs.fasta```

     fasta file containing the representative sequences of the post-CLARC accessory/core COG definitions. The entries in this file will correspond to the COGs in the                                  ```clarc_population_acc_presence_absence.csv``` and ```clarc_population_core_presence_absence.csv``` output files. The COG name will be the fasta entry identifier, with new clusters            having the appropiate cluster name. For clusters, the longest sequence between the individual representative sequences of each COG is chosen as the condensed cluster representative              sequence.
     
   - ```accessory_cluster_cogs.txt``` and ```core_cluster_cogs.txt```

     txt files that include a list of the unique original COGs present in the accessory and core clusters identified by CLARC. No information on the actual clusters (which COGs are connected to      which) is included here, just the list of the original COGs that belong to any cluster.

   - ```accessory_cluster_summary.csv``` and ```core_cluster_summary.csv```

   csv file specifying the compositions of each accessory and core cluster identified by CLARC. The number of columns will depend on the number of members in the largest cluster. The original      frequency in the population of each COG member is also included (freq_cogX) as well as the new frequency of the re-defined cluster (freq_sum).
   
Intermediate output files are:

In the filtering step (before CLARC re-definition), 5 files are created:

   - ```population_accessory_presence_absence.csv```

      csv file containing the presence absence matrix for the _accessory_ genes within the samples specified in the ```needed_sample_names.txt``` input file. So, only COGs present in between          5-95% of those samples will be included here (if run on the default parameters). This is just a filtered output from the original pangenome analysis, before running the CLARC clustering         algorithm.
     
   - ```accessory_rep_seqs.fasta```

      fasta file containing the representative DNA sequences of all _accessory_ COGs identified in the given population of interest. The entries in this file will correspond to the COGs in the        ```population_accessory_presence_absence.csv``` output file. The COG name (as given by the pangenome tool) will be the fasta entry identifier. 

   - ```accessory_rep_protein_seqs.fasta```

      fasta file containing the representative protein sequences of all accessory COGs identified in the given population of interest. It is just the sequences in the                                  ```accessory_rep_seqs.fasta``` file translated from DNA to protein. The pipeline will internally use these protein sequences to functionally annotate the accessory COGs using the EggNOG         database.

   - ```population_core_presence_absence.csv```

      csv file containing the presence absence matrix for the _core_ genes within the samples specified in the ```needed_sample_names.txt``` input file. So, only COGs present in over 95% of           those samples will be included here (if run on default parameters). This is just a filtered output from the original pangenome analysis, before running the CLARC clustering algorithm.           Essentially it is the ```population_accessory_presence_absence.csv``` output, but for core genes instead of accessory.

   - ```core_rep_seqs.fasta```

      fasta file containing the representative DNA sequences of all _core_ COGs identified in the given population of interest. Likewise, this file is equivalent to the                                ```accessory_rep_seqs.fasta``` output file, but for core genes instead of accessory.

### linkage folder

   - ```acc_pXX_matrix.csv```

     csv files that are pairwise matrices with counts of genomes that show the different states of co-occurrence between every possible pair of accessory COGs. The possible states are: genomes       where both COGs are observed together (P11), genomes where neither COG is observed (P00), or genomes where only one of the COGs is observed (P10 and P01). Thus, 4 of these matrices are          generated, one per co-occurrence state.

   - ```acc_linkage_co-occur.csv```

     csv file containing a pairwise metric of linkage. It is internally calculated using the counts in the previous ```acc_pXX_matrix.csv``` matrices, with the following equation:

     <img width="206" alt="image" src="https://github.com/IndraGonz/CLARC/assets/70214497/49f8b26f-f9b9-4463-b12f-92c7d1895f58">

     This provides a normalized linkage metric where:
     
        - L > 0 implies a positive correlation between the genes in the pair (co-occur often)
        - L < 0 implies a negative correlation between the genes in the pair (exclude each other)
        - L = 0 implies no correlation between the genes in the pair

### acc_blastn folder

Various files are created here, since a BLAST database is built from the accessory cog sequence fasta file. However the only relevant file is the one containing the results for the all vs. all accessory COG BLAST:

   - ```blastn_acccogs_allvall.tsv```

     tsv file of the BLAST results with 12 columns, with the following headers:

      ```python
     ['query_seq_ID', 'subject_seq_ID', 'percentage_indentical_matches', 'align_length','num_mismatches','gap_open','align_start_query','align_end_query','align_start_subject', 'align_end_subject','e-value','bit_score']
      ```
      
### eggnog folder

   - ```acc_cog_eggnog_annotations.csv```

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

   - ```eggnog_group_cog_names.csv```

     csv containing a summary of the COGs found in each category. The column headers are the categories, and each column has the list of COGs identified in that category. 

These functional annotation files (and all other intermediate files) are created for the accessory genes filtered from the original pangenome analysis results. However, one of the conditions that CLARC uses to create the same gene clusters is that the members of the cluster must have the same functional group. So, using these annotations the user can also obtain the functional annotated for the new CLARC accessory gene definitions.

## Raising issues

If you have trouble with any of the steps shown here, please let me know! You can email me at igonzalezojeda@g.harvard.edu or submit an issue through GitHub. Also feel free to reach out with any general questions/suggestions!

## Citing

 Linkage-based ortholog refinement in bacterial pangenomes with CLARC
Indra Gonzalez Ojeda, Samantha G Palace, Pamela P Martinez, Taj Azarian, Lindsay R Grant, Laura Hammitt, Bill Hanage, Marc Lipsitch
bioRxiv 2024.12.18.629228; doi: https://doi.org/10.1101/2024.12.18.629228 

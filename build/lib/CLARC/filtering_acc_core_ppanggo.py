#!/usr/bin/env python
# coding: utf-8

#This code has a function that takes the 'gene_presence_absence.Rtab' output from PPanGGolin and filters to obtain a presence absence matrix of only the accessory genes of a population of samples, where the columns are the accessory genes and the indeces are the accession numbers of the samples selected for the analysis. It reads a text file containing the sample names that are to be selected. Accessory genes are defined as those in 5-95% of the selected samples.

#Additionally, from the 'all_nucleotide_families.fasta' file it creates a fasta file with the representative sequences of only the accessory genes.

#There is also a separate function that does this but for the core genes, which are defined as those in >95% of the samples in the population.

### Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
from Bio import SeqIO
import os

### Define functions

## Accessory gene filtering

def get_pop_acc_pres_abs_ppanggo(data_path, out_path, acc_upper, acc_lower):

    # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)
    pres_abs_path = data_path+"/gene_presence_absence.Rtab"

    # Import Ppanggo output
    igopan_all_ppanggo = pd.read_csv(pres_abs_path, sep='\t', comment='#')
    #igopan_all_ppanggo  = igopan_all_ppanggo .T

    # Bakta often has spaces and special characters in its gene names, which can cause problems in the downstream analyses. So here we can turn them into underscores. I'm still deciding if this update is necessary. I think as long as the fasta files pick up the whole name, it should be fine.
    # Update: this was absolutely necessary
    igopan_all_ppanggo["Gene"] = igopan_all_ppanggo["Gene"].str.replace("[ ,\'\"]", "_", regex=True)

    # Get list of PPanGGolin output names in a list
    panppanggo_ids_list =  list(igopan_all_ppanggo["Gene"])

    # Now make gene names the indeces
    #roary_onefilt.set_index('Gene', inplace=True)
    igopan_all_ppanggo.set_index('Gene', inplace=True)

    # Switch rows to columns for ease
    ppanggo_genefreq_matrix = igopan_all_ppanggo.transpose()
    ppanggo_genefreq_matrix.index.name='Accession'

    default_samples = list(igopan_all_ppanggo.columns)
    # Make it so that if no needed_sample_names.txt file is given, it is created based on all isolate names found in the pangenome analysis
    sample_needed_path = data_path+'/needed_sample_names.txt'

    # Check if the file exists
    if not os.path.isfile(sample_needed_path):
        with open(sample_needed_path, 'w') as f:
            for sample in default_samples:
                f.write(sample + '\n')
        print(f"No population sample list specified by user; CLARC used all isolates found in input data to generate the file: {sample_needed_path}")
    else:
        print(f"Population samples specified by user: {sample_needed_path}")

    # Get list of sample names in the subpopulation to analyse
    with open(sample_needed_path, 'r') as file:
        acc_needed_list = file.read().splitlines()

    # Get only isolates in the list
    genefreq_mat_filt = ppanggo_genefreq_matrix[ppanggo_genefreq_matrix.index.isin(acc_needed_list)]

    # Write function that loops through each gene column and returns the frequency in the dataframe
    def get_freq(dataframe):

        cog_list = []
        freq_list = []

        for element in dataframe.columns:

            freq_dataframe = pd.DataFrame()
            cog_freq = (dataframe[element].sum())/(dataframe.shape[0])

            cog_list.append(element)
            freq_list.append(cog_freq)

        freq_dataframe['COG_name'] = cog_list
        freq_dataframe['freq'] = freq_list

        return(freq_dataframe)


    # Get list of genes and their frequency for desired population

    #Navajo all
    freq_cog_navajo = get_freq(genefreq_mat_filt)

    # Make sure the frequencies are numeric
    freq_cog_navajo['freq'] = pd.to_numeric(freq_cog_navajo['freq'], errors='coerce')

    # Now filter to only keep cogs with frequency between 5 and 95% within the given sample subset
    acc_cog_navajo = freq_cog_navajo[(freq_cog_navajo['freq'] > acc_lower) & (freq_cog_navajo['freq'] < acc_upper)]

    # Get the number of accessory cogs per dataset
    acccog_num_navajo = acc_cog_navajo.shape[0]

    # Get accessory cog names as lists
    acccog_name_navajo_list = list(acc_cog_navajo["COG_name"])

    # Only get presence absence matrix for accessory genes
    genefreq_meta_filt_acc = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(acccog_name_navajo_list)]]

    ## Export results to output folder
    csv_path = out_path+f'/population_accessory_presence_absence.csv'

    genefreq_meta_filt_acc.to_csv(csv_path, index=True)

    ### Now get fasta file with only accessory genes

    # Define the input and output file paths
    seq_path = data_path+"/all_nucleotide_families.fasta"
    acc_fasta_path = out_path+f"/accessory_rep_seqs.fasta"

    # Open the input and output files
    with open(seq_path, "r") as infile, open(acc_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id = record.description
            cog_id_fixed = cog_id.replace(' ', '_').replace(',', '_').replace("'", '_').replace('"', '_')
            if cog_id_fixed in acccog_name_navajo_list:
                # Write the record to the output file
                record.id = cog_id_fixed
                record.description = ''
                SeqIO.write(record, outfile, "fasta")

## Core gene filtering

def get_pop_core_pres_abs_ppanggo(data_path, out_path, core_lower):

    # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)
    pres_abs_path = data_path+"/gene_presence_absence.Rtab"
    sample_needed_path = data_path+'/needed_sample_names.txt'

    # Import Ppanggo output
    igopan_all_ppanggo = pd.read_csv(pres_abs_path, sep='\t', comment='#')
    #igopan_all_ppanggo  = igopan_all_ppanggo .T

    # Bakta often has spaces and special characters in its gene names, which can cause problems in the downstream analyses. So here we can turn them into underscores. I'm still deciding if this update is necessary. I think as long as the fasta files pick up the whole name, it should be fine.
    # Update: this was absolutely necessary
    igopan_all_ppanggo["Gene"] = igopan_all_ppanggo["Gene"].str.replace("[ ,\'\"]", "_", regex=True)

    # Get list of PPanGGolin output names in a list
    panppanggo_ids_list =  list(igopan_all_ppanggo["Gene"])

    # Now make gene names the indeces
    #roary_onefilt.set_index('Gene', inplace=True)
    igopan_all_ppanggo.set_index('Gene', inplace=True)

    # Switch rows to columns for ease
    ppanggo_genefreq_matrix = igopan_all_ppanggo.transpose()
    ppanggo_genefreq_matrix.index.name='Accession'

    # Get list of sample names in the subpopulation to analyse
    with open(sample_needed_path, 'r') as file:
        acc_needed_list = file.read().splitlines()

    # Get only isolates in the list
    genefreq_mat_filt = ppanggo_genefreq_matrix[ppanggo_genefreq_matrix.index.isin(acc_needed_list)]

    # Write function that loops through each gene column and returns the frequency in the dataframe
    def get_freq(dataframe):

        cog_list = []
        freq_list = []

        for element in dataframe.columns:

            freq_dataframe = pd.DataFrame()
            cog_freq = (dataframe[element].sum())/(dataframe.shape[0])

            cog_list.append(element)
            freq_list.append(cog_freq)

        freq_dataframe['COG_name'] = cog_list
        freq_dataframe['freq'] = freq_list

        return(freq_dataframe)


    # Get list of genes and their frequency for desired population

    #Navajo all
    freq_cog_navajo = get_freq(genefreq_mat_filt)

    # Make sure the frequencies are numeric
    freq_cog_navajo['freq'] = pd.to_numeric(freq_cog_navajo['freq'], errors='coerce')

    # Now filter to only keep cogs with frequency over a given frequency (default 95%) within the given sample subset
    core_cog_navajo = freq_cog_navajo[(freq_cog_navajo['freq'] >= core_lower)]

    # Get the number of core cogs per dataset
    corecog_num_navajo = core_cog_navajo.shape[0]

    # Get core cog names as lists
    corecog_name_navajo_list = list(core_cog_navajo["COG_name"])

    # Only get presence absence matrix for core genes
    genefreq_meta_filt_core = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(corecog_name_navajo_list)]]

    ## Export results to output folder
    csv_path = out_path+f'/population_core_presence_absence.csv'

    genefreq_meta_filt_core.to_csv(csv_path, index=True)

    ### Now get fasta file with only core genes

    # Define the input and output file paths
    seq_path = data_path+"/all_nucleotide_families.fasta"
    core_fasta_path = out_path+f"/core_rep_seqs.fasta"

    # Open the input and output files
    with open(seq_path, "r") as infile, open(core_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id = record.description
            cog_id_fixed = cog_id.replace(' ', '_').replace(',', '_').replace("'", '_').replace('"', '_')
            if cog_id_fixed in corecog_name_navajo_list:
                # Write the record to the output file
                record.id = cog_id_fixed
                record.description = ''
                SeqIO.write(record, outfile, "fasta")

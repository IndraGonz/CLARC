#!/usr/bin/env python
# coding: utf-8

### Filtering pangenome presence absence matrix to get accessory and core genes

# This code has a function that takes the 'gene_presence_absence.csv' output from Roary and filters to obtain a presence absence matrix of only the accessory genes of a population of samples, where the columns are the accessory genes and the indeces are the accession numbers of the samples selected for the analysis. It reads a text file containing the sample names that are to be selected. Accessory genes as default defined as those in 5-95% of the selected samples, but bound can be specified otherwise.
#
# Additionally, from the 'pan_genome_reference.fa' file it creates a fasta file with the representative sequences of only the accessory genes.
#
# There is also a separate function that does this but for the core genes, which in the default are defined as those in >95% of the samples in the population, but a different lower bound can be specified.

### Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
from Bio import SeqIO

### Define functions

# Function to filter for accessory cogs

def get_pop_acc_pres_abs(data_path, out_path, acc_upper, acc_lower):

    # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)

    pres_abs_path = data_path+"/gene_presence_absence.csv"
    sample_needed_path = data_path+'/needed_sample_names.txt'

    # Import Roary output
    igopan_all_roary = pd.read_csv(pres_abs_path, low_memory=False)

    # Bakta often has spaces and special characters in its gene names, which can cause problems in the downstream analyses. So here we can turn them into underscores. I'm still deciding if this update is necessary. I think as long as the fasta files pick up the whole name, it should be fine.
    # Update: this was absolutely necessary
    igopan_all_roary["Gene"] = igopan_all_roary["Gene"].str.replace("[ ,\'\"]", "_", regex=True)

    # Get list of Roary output names in a list
    panroary_ids_list =  list(igopan_all_roary["Gene"])

    # I have eliminated this filtering step, since we were using it for S. pneumoniae, but it might not hold for all species and it is unclear how necessary it is
    # First round of filtering by fragment length >150bp OR general isolate frequency >10%
    #tenp = (igopan_all_roary.shape[1]-14)/10 # Counting number of isolates (first 14 columns are metadata)
    #roary_onefilt = igopan_all_roary[(igopan_all_roary['Avg group size nuc'] >= 150) | (igopan_all_roary['No. isolates'] >= tenp)]

    # Now make gene names the indeces
    #roary_onefilt.set_index('Gene', inplace=True)
    igopan_all_roary.set_index('Gene', inplace=True)

    # Drop all columns that are not an isolate
    #roary_isol = roary_onefilt.iloc[:,13:]
    roary_isol = igopan_all_roary.iloc[:,13:]

    # Now replace NaN values for 0 and any other value for 1
    roary_isol[~roary_isol.isnull()] = 1
    roary_isol[roary_isol.isnull()] = 0

    # Switch rows to columns for ease
    roary_genefreq_matrix = roary_isol.transpose()
    roary_genefreq_matrix.index.name='Accession'

    # Get list of sample names in the subpopulation to analyse
    with open(sample_needed_path, 'r') as file:
        acc_needed_list = file.read().splitlines()

    # Get only isolates in the list
    genefreq_mat_filt = roary_genefreq_matrix[roary_genefreq_matrix.index.isin(acc_needed_list)]

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
    acc_cog_navajo = freq_cog_navajo[(freq_cog_navajo['freq'] >= acc_lower) & (freq_cog_navajo['freq'] <= acc_upper)]

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
    seq_path = data_path+"/pan_genome_reference.fa"
    acc_fasta_path = out_path+f"/accessory_rep_seqs.fasta"

    # Open the input and output files
    with open(seq_path, "r") as infile, open(acc_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id = record.description.split(' ', 1)[1]
            cog_id_fixed = cog_id.replace(' ', '_').replace(',', '_').replace("'", '_').replace('"', '_')
            if cog_id_fixed in acccog_name_navajo_list:
                # Write the record to the output file
                record.id = cog_id_fixed
                record.description = ''
                SeqIO.write(record, outfile, "fasta")


### Function to filter for core cogs

def get_pop_core_pres_abs(data_path, out_path, core_lower):

    # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)

    pres_abs_path = data_path+"/gene_presence_absence.csv"
    sample_needed_path = data_path+'/needed_sample_names.txt'

    # Import Roary output
    igopan_all_roary = pd.read_csv(pres_abs_path, low_memory=False)

    # Bakta often has spaces and special characters in its gene names, which can cause problems in the downstream analyses. So here we turn them into underscores.
    igopan_all_roary["Gene"] = igopan_all_roary["Gene"].str.replace("[ ,\'\"]", "_", regex=True)

    # Get list of Roary output names in a list
    panroary_ids_list =  list(igopan_all_roary["Gene"])

    # I have eliminated this filtering step, since we were using it for S. pneumoniae, but it might not hold for all species and it is unclear how necessary it is
    # First round of filtering by fragment length >150bp OR general isolate frequency >10%
    #tenp = (igopan_all_roary.shape[1]-14)/10 # Counting number of isolates (first 14 columns are metadata)
    #roary_onefilt = igopan_all_roary[(igopan_all_roary['Avg group size nuc'] >= 150) | (igopan_all_roary['No. isolates'] >= tenp)]

    # Now make gene names the indeces
    #roary_onefilt.set_index('Gene', inplace=True)
    igopan_all_roary.set_index('Gene', inplace=True)

    # Drop all columns that are not an isolate
    #roary_isol = roary_onefilt.iloc[:,13:]
    roary_isol = igopan_all_roary.iloc[:,13:]

    # Now replace NaN values for 0 and any other value for 1
    roary_isol[~roary_isol.isnull()] = 1
    roary_isol[roary_isol.isnull()] = 0

    # Switch rows to columns for ease
    roary_genefreq_matrix = roary_isol.transpose()
    roary_genefreq_matrix.index.name='Accession'

    # Get list of sample names in the subpopulation to analyse
    with open(sample_needed_path, 'r') as file:
        acc_needed_list = file.read().splitlines()

    # Get only isolates in the list
    genefreq_mat_filt = roary_genefreq_matrix[roary_genefreq_matrix.index.isin(acc_needed_list)]

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
    core_cog_navajo = freq_cog_navajo[(freq_cog_navajo['freq'] > core_lower)]

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
    seq_path = data_path+"/pan_genome_reference.fa"
    core_fasta_path = out_path+f"/core_rep_seqs.fasta"

    # Open the input and output files
    with open(seq_path, "r") as infile, open(core_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id = record.description.split(' ', 1)[1]
            cog_id_fixed = cog_id.replace(' ', '_').replace(',', '_').replace("'", '_').replace('"', '_')
            if cog_id_fixed in corecog_name_navajo_list:
                # Write the record to the output file
                record.id = cog_id_fixed
                record.description = ''
                SeqIO.write(record, outfile, "fasta")

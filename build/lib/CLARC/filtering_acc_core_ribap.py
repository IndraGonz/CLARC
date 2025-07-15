#!/usr/bin/env python
# coding: utf-8

# # Filtering pangenome presence absence matrix to get accessory and core genes - RIBAP output

# **Here I will revise CLARC's latest version (1.2.1) to implement additional features, including possibly adding RIBAP output as an option**

# This code has a function that takes the 'holy_python_ribap_95.csv' output from RIBAP and filters to obtain a presence absence matrix of only the accessory genes of a population of samples, where the columns are the accessory genes and the indeces are the accession numbers of the samples selected for the analysis. It reads a text file containing the sample names that are to be selected. In default parameters, accessory genes are defined as those in 5-95% of the selected samples.
#
# Additionally, from the 'pan_genome_reference.fa' file in the '/02-roary/95' folder within the RIBAP results folder, it creates both a global 'pan_genome_reference_ribap.fa' with all RIBAP clusters and a fasta file with the representative sequences of only the accessory genes.
#
# There is also a separate function that does this but for the core genes, which are defined as those in >95% of the samples in the population.

# ### Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import time
import os

# ### Define functions

def get_pop_acc_pres_abs_ribap(data_path, out_path, acc_upper, acc_lower):

    # Import RIBAP output
    pres_abs_path = data_path+'/05-combine/holy_python_ribap_95.csv'

    igopan_all_ribap = pd.read_csv(pres_abs_path, sep="\t")
    igopan_all_ribap["Cluster_ID"] = igopan_all_ribap["Cluster_ID"].str.replace("[ ,\'\"]", "_", regex=True)
    panribap_ids_list =  list(igopan_all_ribap["Cluster_ID"])

    # Now make gene names the indeces
    igopan_all_ribap.set_index('Cluster_ID', inplace=True)

    # Drop all columns that are not an isolate
    ribap_isol = igopan_all_ribap.iloc[:,2:]

    ribap_isol_tags = ribap_isol.copy()

    # Now replace NaN values for 0 and any other value for 1
    ribap_isol[~ribap_isol.isnull()] = 1
    ribap_isol[ribap_isol.isnull()] = 0

    # Switch rows to columns for ease
    ribap_genefreq_matrix = ribap_isol.transpose()
    ribap_genefreq_matrix.index.name='Accession'

    default_samples = list(ribap_genefreq_matrix.index)

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
    genefreq_mat_filt = ribap_genefreq_matrix[ribap_genefreq_matrix.index.isin(acc_needed_list)]

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
    acc_cog_navajo = freq_cog_navajo[(freq_cog_navajo['freq'] >= acc_lower) & (freq_cog_navajo['freq'] < acc_upper)]

    # Get the number of accessory cogs per dataset
    acccog_num_navajo = acc_cog_navajo.shape[0]

    # Get accessory cog names as lists
    acccog_name_navajo_list = list(acc_cog_navajo["COG_name"])

    #### Next step: Generate pan_genome_reference_ribap.fa

    ribap_roary_fasta_path = data_path+'/02-roary/95/pan_genome_reference.fa'
    ribap_out_fasta_path = data_path+'/pan_genome_reference_ribap.fa'
    fasta_dict = SeqIO.to_dict(SeqIO.parse(ribap_roary_fasta_path, "fasta"))

    output_records = []

    for group_id, row in ribap_isol_tags.iterrows():
        sequences = []

        for tag in row.dropna():
            if tag != 'NA' and tag in fasta_dict:
                sequences.append(fasta_dict[tag])

        if sequences:
            # Select the longest sequence
            longest_seq = max(sequences, key=lambda x: len(x.seq))
            new_record = SeqRecord(
                longest_seq.seq,
                id=group_id,
                description=""
            )
            output_records.append(new_record)
        else:
            print(f"Warning: No sequence found for group '{group_id}'")

    SeqIO.write(output_records, ribap_out_fasta_path, "fasta")

    #### Continue filtering by creating the presence absence matrix only for the accessory COGs

    # Only get presence absence matrix for accessory genes
    genefreq_meta_filt_acc = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(acccog_name_navajo_list)]]

    ## Export results to output folder
    csv_path = out_path+f'/population_accessory_presence_absence.csv'

    genefreq_meta_filt_acc.to_csv(csv_path, index=True)

    ### Now get fasta file with only accessory genes

    # Define the input and output file paths
    seq_path = data_path+"/pan_genome_reference_ribap.fa"
    acc_fasta_path = out_path+f"/accessory_rep_seqs.fasta"

    # Open the input and output files
    with open(seq_path, "r") as infile, open(acc_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id_fixed = record.id
            if cog_id_fixed in acccog_name_navajo_list:
                # Write the record to the output file
                record.id = cog_id_fixed
                record.description = ''
                SeqIO.write(record, outfile, "fasta")


def get_pop_core_pres_abs_ribap(data_path, out_path, core_lower):

    # Import RIBAP output
    pres_abs_path = data_path+'/05-combine/holy_python_ribap_95.csv'

    igopan_all_ribap = pd.read_csv(pres_abs_path, sep="\t")
    igopan_all_ribap["Cluster_ID"] = igopan_all_ribap["Cluster_ID"].str.replace("[ ,\'\"]", "_", regex=True)
    panribap_ids_list =  list(igopan_all_ribap["Cluster_ID"])

    # Now make gene names the indeces
    igopan_all_ribap.set_index('Cluster_ID', inplace=True)

    # Drop all columns that are not an isolate
    ribap_isol = igopan_all_ribap.iloc[:,2:]

    ribap_isol_tags = ribap_isol.copy()

    # Now replace NaN values for 0 and any other value for 1
    ribap_isol[~ribap_isol.isnull()] = 1
    ribap_isol[ribap_isol.isnull()] = 0

    # Switch rows to columns for ease
    ribap_genefreq_matrix = ribap_isol.transpose()
    ribap_genefreq_matrix.index.name='Accession'

    default_samples = list(ribap_genefreq_matrix.index)

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
    genefreq_mat_filt = ribap_genefreq_matrix[ribap_genefreq_matrix.index.isin(acc_needed_list)]

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
    core_cog_navajo = freq_cog_navajo[(freq_cog_navajo['freq'] >= core_lower)]

    # Get the number of core cogs per dataset
    corecog_num_navajo = core_cog_navajo.shape[0]

    # Get core cog names as lists
    corecog_name_navajo_list = list(core_cog_navajo["COG_name"])

    #### Continue filtering by creating the presence absence matrix only for the core COGs

    # Only get presence absence matrix for core genes
    genefreq_meta_filt_core = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(corecog_name_navajo_list)]]

    ## Export results to output folder
    csv_path = out_path+f'/population_core_presence_absence.csv'

    genefreq_meta_filt_core.to_csv(csv_path, index=True)

    ### Now get fasta file with only core genes

    # Define the input and output file paths
    seq_path = data_path+"/pan_genome_reference_ribap.fa"
    core_fasta_path = out_path+f"/core_rep_seqs.fasta"

    # Open the input and output files
    with open(seq_path, "r") as infile, open(core_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id_fixed = record.id
            if cog_id_fixed in corecog_name_navajo_list:
                # Write the record to the output file
                record.id = cog_id_fixed
                record.description = ''
                SeqIO.write(record, outfile, "fasta")

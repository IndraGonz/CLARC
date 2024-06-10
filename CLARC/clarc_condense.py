#!/usr/bin/env python
# coding: utf-8

# # Condensing COG definitions with CLARC method

# Here we use linkage, functional annotations and sequence identity to identify clusters of COGs that are the same gene.

import pandas as pd
import numpy as np
import scipy as scipy
import itertools
import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt
import fastcluster
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


# # Import linkage results

def clarc_cleaning(out_path):

    # Get list of COG pairs that do not co-occur within the dataset

    # Import counts for COG pairs that co-occur
    linkage_p11_path = out_path+"/linkage/acc_p11_matrix.csv"
    linkage_p11_result = pd.read_csv(linkage_p11_path ,low_memory=False, index_col=0)

    # Import log co-occurrence matrix
    linkage_logcoccur_path = out_path+"/linkage/acc_linkage_co-occur.csv"
    linkage_logcoccur_result = pd.read_csv(linkage_logcoccur_path, low_memory=False, index_col=0)

    # Check COG pairs that never co-occur

    p11_coglist = linkage_p11_result.melt(ignore_index = False)
    p11_coglist = p11_coglist.reset_index()
    p11_coglist.columns = ['COG1','COG2','p11']

    # Turning the log co-occur matrix into list
    linkage_logcoccur_list = linkage_logcoccur_result.melt(ignore_index = False)
    linkage_logcoccur_list = linkage_logcoccur_list.reset_index()
    linkage_logcoccur_list.columns = ['COG1','COG2','log_cooccur']

    p11_coglist_zeros = p11_coglist[p11_coglist["p11"] == 0]

    # Now we have to drop the duplicates from these dataframes

    # p11
    df_p11 = pd.DataFrame(np.sort(p11_coglist_zeros[['COG1','COG2']].values,1), columns=['COG1', 'COG2']).drop_duplicates()
    p11_coglist_zeros_unique = p11_coglist_zeros.merge(df_p11, on=['COG1', 'COG2'])

    # Calculate COG frequencies in the population

    # Import presence absence matrix
    pres_abs_path = out_path+'/population_accessory_presence_absence.csv'
    pan_acc = pd.read_csv(pres_abs_path, index_col=0, low_memory=False)

    # Calculate freqs of each COG in subpopulation dataset
    cog_freq = pan_acc.mean()
    cog_freq = cog_freq.transpose()
    cog_freq_df = cog_freq.to_frame()
    cog_freq_df = cog_freq_df.reset_index()
    cog_freq_df.columns = ['COG','freq']
    # Generate colums to be merged with cog pair dataframe
    cog_freq_df_cog1 = cog_freq_df.copy()
    cog_freq_df_cog1.columns = ['COG1','freq_cog1']
    cog_freq_df_cog2 = cog_freq_df.copy()
    cog_freq_df_cog2.columns = ['COG2','freq_cog2']

    # Import eggnog functional annotations

    funct_path = out_path+'/eggnog/eggnog_group_cog_names.csv'
    func_subset_names = pd.read_csv(funct_path)

    # eggnog
    eggnog_func_subsets = func_subset_names.drop(columns=['all_acc_COGs'])
    eggnog_cog_category_df = eggnog_func_subsets.melt()
    eggnog_cog_category_df = eggnog_cog_category_df.dropna(axis='index')
    eggnog_cog_category_df = eggnog_cog_category_df.rename(columns={"value": "cog", "variable": "functional_group"})
    column_order = ['cog'] + [col for col in eggnog_cog_category_df.columns if col != 'cog']
    eggfunc_cog = eggnog_cog_category_df[column_order]

    # Generate colums to be merged with cog pair dataframe
    eggfunc_cog1 = eggfunc_cog.copy()
    eggfunc_cog1.columns = ['COG1','funct_group_cog1']
    eggfunc_cog2 = eggfunc_cog.copy()
    eggfunc_cog2.columns = ['COG2','funct_group_cog2']

    # Merge
    p11_coglist_freq = p11_coglist_zeros_unique.copy()
    p11_coglist_freq = p11_coglist_freq.merge(cog_freq_df_cog1,on='COG1')
    p11_coglist_freq = p11_coglist_freq.merge(cog_freq_df_cog2,on='COG2')
    p11_coglist_freq = p11_coglist_freq.merge(eggfunc_cog1,on='COG1')
    p11_coglist_freq = p11_coglist_freq.merge(eggfunc_cog2,on='COG2')

    # Calculate samples nedeed to see one co-occurrence by chance and other things
    p11_coglist_freq['P(12)'] = p11_coglist_freq['freq_cog1']*p11_coglist_freq['freq_cog2']
    p11_coglist_freq['freq_cog1+freq_cog2'] = p11_coglist_freq['freq_cog1']+p11_coglist_freq['freq_cog2']
    p11_coglist_freq['freq_cog1-freq_cog2'] = p11_coglist_freq['freq_cog1']-p11_coglist_freq['freq_cog2']

    # Add flag to check if the eggnog func group is the same
    p11_coglist_freq.loc[p11_coglist_freq['funct_group_cog1']==p11_coglist_freq['funct_group_cog2'], 'eggnog_same_flag'] = 1

    # Adding blastn results for the accessory COGs against each other
    blast_path = out_path+'/acc_blastn/blastn_acccogs_allvall.tsv'
    blastn_acccogs = pd.read_table(blast_path, header=None)

    # TSV output doesn't have column labels, so here I add them
    column_labels = ['query_seq_ID', 'subject_seq_ID', 'percentage_identical_matches', 'align_length','num_mismatches','gap_open','align_start_query','align_end_query','align_start_subject', 'align_end_subject','e-value','bit_score']
    blastn_acccogs.columns = column_labels

    # Add empty column to dataframe with cog pairs that exclude each other
    p11_coglist_freq['blastn_iden_per'] = ""

    for row in p11_coglist_freq.iterrows():

        # Get COG names for row
        cog1 = row[1].COG1
        cog2 = row[1].COG2

        # Search in blast results dataframe
        blast_cog1cog2 = blastn_acccogs.query(f"query_seq_ID == '{cog1}' and subject_seq_ID == '{cog2}'")


        # Check if it does have a hit
        if blast_cog1cog2.empty:
            p11_coglist_freq.at[row[0], 'blastn_iden_per'] = np.nan
        else:
            p11_coglist_freq.at[row[0], 'blastn_iden_per'] = blast_cog1cog2['percentage_identical_matches'].values[0]

    # Same eggnog functional group
    samecog_condense = p11_coglist_freq[p11_coglist_freq['eggnog_same_flag'] == 1]

    # High blastn sequence identity
    samecog_condense = samecog_condense[samecog_condense['blastn_iden_per'] >= 95]

    # Get list of unique COGs

    p11_zero_cog1_list = list(samecog_condense['COG1'])
    p11_zero_cog2_list = list(samecog_condense['COG2'])

    samegene_cogs = p11_zero_cog1_list + p11_zero_cog2_list
    samegene_cogs_list = set(samegene_cogs)
    samegene_cogs_list = list(samegene_cogs_list)

    # Loop through each COG that appears in the dataframe where COGs are in >1 pair

    samegene_connected_cogs = pd.DataFrame(columns=['connected_COGs','group_number','group_length'])
    group = 0

    for cog in samegene_cogs_list:

        # Get COGs with which that COG has a relationship
        con = samecog_condense.query(f"COG1 == '{cog}' or COG2 == '{cog}'") # Connections of the COG
        x = sorted(pd.concat([con['COG1'],con['COG2']]).unique())
        xl = sorted(pd.concat([con['COG1'],con['COG2']]).unique())
        xl.remove(cog)

        # Here I do a loop where if the COGs that have relationships to the main COG ONLY have relationships to each other, then it is added to the dataframe

        flag = 0

        for i in xl:

            coni = samecog_condense.query(f"COG1 == '{i}' or COG2 == '{i}'") # Connections of the COG
            y = sorted(pd.concat([coni['COG1'],coni['COG2']]).unique())

            if y == x:

                flag = flag + 1

        if flag == len(xl):

            group = group + 1

            samegene_connected_cogs.loc[group-1,'connected_COGs'] = x
            samegene_connected_cogs.loc[group-1,'group_number'] = group
            samegene_connected_cogs.loc[group-1,'group_length'] = len(x)

    samegene_connected_cogs['connected_COGs'] = samegene_connected_cogs['connected_COGs'].apply(tuple)
    samegene_connected_cogs = samegene_connected_cogs.drop_duplicates(subset='connected_COGs', keep="first")
    samegene_connected_cogs = samegene_connected_cogs.sort_values(by=['group_length'], ascending=False)
    samegene_connected_cogs = samegene_connected_cogs.reset_index(drop=True)

    def summarize_clusters(samegene_connected_cogs, cog_freq_df, max_clusters):
        # Create dataframe to summarize the clusters
        columns = [f'COG{i}' for i in range(1, max_clusters+1)]
        samegene_connected_cogs_exp = pd.DataFrame(samegene_connected_cogs.iloc[:, 0].to_list(), columns=columns)

        # Set empty columns where frequency values will go
        for i in range(1, max_clusters+1):
            samegene_connected_cogs_exp[f'freq_cog{i}'] = ""

        # Get frequency values
        for i in range(1, max_clusters+1):
            cog_list = samegene_connected_cogs_exp[f'COG{i}']
            for j, cog in enumerate(cog_list):
                if cog is None:
                    samegene_connected_cogs_exp.at[j, f'freq_cog{i}'] = None
                else:
                    x = cog_freq_df.query(f"COG == '{cog}'")
                    x = x.reset_index()
                    freq = x['freq'][0] if not x.empty else None
                    samegene_connected_cogs_exp.at[j, f'freq_cog{i}'] = freq

        return samegene_connected_cogs_exp

    max_clusters = samegene_connected_cogs['group_length'][0]
    samegene_connected_cogs_exp = summarize_clusters(samegene_connected_cogs, cog_freq_df, max_clusters)

    # Get the list of column names
    frequency_columns = [f"freq_cog{i}" for i in range(1, max_clusters+1)]

    # Add freq sum
    samegene_connected_cogs_exp["freq_sum"] = samegene_connected_cogs_exp[frequency_columns].sum(axis=1)

    # Now get clusters that seem to be core genes and clusters that seem to be accessory genes

    # They are core if their joint frequency is >= 0.95
    conn_cluster_cogs_core = samegene_connected_cogs_exp[samegene_connected_cogs_exp['freq_sum'] >= 0.95]

    # They are accessory is their joint frequency is < 0.95
    conn_cluster_cogs_acc = samegene_connected_cogs_exp[samegene_connected_cogs_exp['freq_sum'] < 0.95]

    cog_columns = [f"COG{i}" for i in range(1, max_clusters + 1)]

    cluster_corecogs_list = []
    for col in cog_columns:
        cog_list_a = list(conn_cluster_cogs_core[col])
        cog_list = list(filter(None, cog_list_a))
        cluster_corecogs_list.extend(cog_list)

    cluster_corecogs_unique = set(cluster_corecogs_list)

    cog_columns = [f"COG{i}" for i in range(1, max_clusters + 1)]

    cluster_acccogs_list = []
    for col in cog_columns:
        cog_list_a = list(conn_cluster_cogs_acc[col])
        cog_list = list(filter(None, cog_list_a))
        cluster_acccogs_list.extend(cog_list)

    cluster_acccogs_unique = set(cluster_acccogs_list)

    clarc_summary_path = out_path+"/clarc_results"

    # Check if the directory already exists (create if it doesn't)
    if not os.path.exists(clarc_summary_path):
        os.makedirs(clarc_summary_path)

    out_text = clarc_summary_path+'/clarc_cluster_summary.txt'
    with open(out_text, "a") as file:
        file.write(f'Core COG clusters: '+f'{len(conn_cluster_cogs_core)}\n')
        file.write(f'Unique COGs in core clusters: '+f'{len(cluster_corecogs_unique)}\n')
        file.write(f'Accessory COG clusters: '+f'{len(conn_cluster_cogs_acc)}\n')
        file.write(f'Unique COGs in accessory clusters: '+f'{len(cluster_acccogs_unique)}\n')

    core_cluster_path = out_path+"/clarc_results/core_cluster_summary.csv"
    conn_cluster_cogs_core.to_csv(core_cluster_path, index=False)

    acc_cluster_path = out_path+"/clarc_results/accessory_cluster_summary.csv"
    conn_cluster_cogs_acc.to_csv(acc_cluster_path, index=False)

    # Eliminate core clusters and join accessory clusters - New presence absence matrix

    conn_cluster_cogs_acc_legend = pd.DataFrame()

    for col in conn_cluster_cogs_acc.columns:
        if 'freq' not in col:
            # Add '-' after each COG value
            conn_cluster_cogs_acc_legend[col] = conn_cluster_cogs_acc[col] + '-'

    # Concatenate COG values with '-' in between
    conn_cluster_cogs_acc_legend['new_cluster_name_samecog'] = conn_cluster_cogs_acc_legend.fillna('').sum(axis=1)

    for col in conn_cluster_cogs_acc.columns:
        if 'freq' not in col:
            conn_cluster_cogs_acc_legend[col] = conn_cluster_cogs_acc[col]

    # Add empty column to dataframe with cog pairs that exclude each other
    pres_abs_newclusters = pd.DataFrame()

    for index, row in conn_cluster_cogs_acc_legend.iterrows():

        # Get COG names for row and new name
        cog_columns = [f'COG{i}' for i in range(1, max_clusters+1)]
        cogs = [row[cog] for cog in cog_columns if pd.notnull(row[cog])]
        new_name = row['new_cluster_name_samecog']

        # Initialize dataframe
        df = pan_acc[cogs].copy()  # Use .copy() to avoid SettingWithCopyWarning

        # Concatenate COGs
        df[new_name] = df.sum(axis=1)
        pres_abs_newclusters = pd.concat([pres_abs_newclusters, df[new_name]], axis=1)

    # Now we want to drop the COGs in the clusters we just condensed and add the presence absence info of the newly defined clusters

    # Drop old redundant COGs from presence absence matrix
    nonredundant_presabs = pan_acc.drop(columns=cluster_acccogs_unique)

    # Now append new groups
    samecog_clustered_presabs = pd.concat([nonredundant_presabs, pres_abs_newclusters], axis=1)

    # Now we drop COGs identified as part of core gene clusters
    samecog_clustered_presabs_nocore = samecog_clustered_presabs.loc[:, ~samecog_clustered_presabs.columns.isin(cluster_corecogs_unique)]

    clarc_presabs_path = out_path+"/clarc_results/clarc_condensed_presence_absence.csv"
    samecog_clustered_presabs_nocore.to_csv(clarc_presabs_path, index=True)

    # Define the file path
    core_list_path = out_path+"/clarc_results/core_cluster_cogs.txt"

    # Open the file in write mode
    with open(core_list_path, "w") as file:
        # Write each element of the list to the file
        for item in cluster_corecogs_unique:
            file.write(item + "\n")  # Add a newline after each item

    # Define the file path
    acc_list_path = out_path+"/clarc_results/accessory_cluster_cogs.txt"

    # Open the file in write mode
    with open(acc_list_path, "w") as file:
        # Write each element of the list to the file
        for item in cluster_acccogs_unique:
            file.write(item + "\n")  # Add a newline after each item

    # Get fasta file with the representative sequences of the new accessory gene definitions

    # Get dataframe with all accessory COGs and their length

    def fasta_coglen_dataframe(fasta_file):

        seq_ids = []
        seq_lengths = []

        for record in SeqIO.parse(fasta_file, "fasta"):

            seq_ids.append(record.id)
            seq_lengths.append(len(record.seq))

        df = pd.DataFrame({"cog_name": seq_ids, "length": seq_lengths})

        return df

    # Get dataframe with length of each accessory cog

    # Provide the path to your fasta file
    original_fasta_path = out_path+'/accessory_rep_seqs.fasta'

    # Parse the fasta file and create the DataFrame
    cog_len = fasta_coglen_dataframe(original_fasta_path)

    conn_cluster_cogs_acc_legend = conn_cluster_cogs_acc_legend.reset_index(drop=True)

    # Identify longest COG in each accessory cluster

    for index, row in conn_cluster_cogs_acc_legend.iterrows():

        longest_cog = ''
        max_len = 0

        for i in range(1, max_clusters):

            cog = row[f'COG{i}']

            if cog:

                len_row = cog_len.query(f"cog_name == '{cog}'")
                leng = len_row['length'].values[0]

                if leng > max_len:

                    longest_cog = cog
                    max_len = leng

            conn_cluster_cogs_acc_legend.at[index, 'longest_cog'] = longest_cog

    # Get important lists
    longest_acc_cogs = list(conn_cluster_cogs_acc_legend['longest_cog'])

    # Get fasta file with the representative sequences for the CLARC COG definitions

    # The COGs identified as part of core clusters will be dropped

    # The COGs identified as part of an accessory cluster, but not the longest cog will be dropped

    # The COGs identified as part of an accessory cluster and representing the longest COG will have their name reannotated with the cluster name

    clarc_fasta_path = out_path+'/clarc_results/clarc_acc_cog_seqs.fasta'

    fin = open(original_fasta_path, 'r')
    fout = open(clarc_fasta_path, 'w')

    for record in SeqIO.parse(fin,'fasta'):

        cog_name = record.id

        if cog_name not in cluster_corecogs_unique:

            if cog_name in cluster_acccogs_unique:

                if cog_name in longest_acc_cogs:

                    x = conn_cluster_cogs_acc_legend.query(f" longest_cog == '{cog_name}'")
                    cluster_name = x['new_cluster_name_samecog'].values[0]
                    record.id = cluster_name
                    SeqIO.write(record, fout, 'fasta')

            else:
                SeqIO.write(record, fout, 'fasta')

    fin.close()
    fout.close()

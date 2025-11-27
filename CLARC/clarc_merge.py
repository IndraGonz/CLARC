#!/usr/bin/env python
# coding: utf-8

# # Merge CLARC cluster results if run on multiple subpopulations of the same pangenome analysis

# In the case of a multipopulation pangenome analysis, CLARC can be run on various subpopulations independently. In this case, the clusters found in each populations can be different.
#
# So, this script serves to unify the results from CLARC run on multiple subpopulations within the same pangenome analysis. This leverages linkage information from distinct populations to further optimize the accessory COG definitions.
#
# The methodology consist of a series of steps to keep clusters with no conflicting information across the populations and eliminate those with conflicting ifnormation (e.g. whehn a COG is part of completely different clusters in two diferent populations).
#
# The steps are as follows:
#
# Step 1: Get CLARC graphs for each distinct population
#
# Step 2: Group clusters in each population based on shared COGs
#
# Step 3: Determine what clusters to keep (keep those with no conflicting information, eliminate those with conflicting information)
#
# Step 4: Get final list of clusters to condense
#
# Step 5: Condense on pangenome output with all populations
#
# The input to this function is the paths of the folders where CLARC was run on each population and the output is the initial presence absence matrix with the CLARC clusters condensed, a summary of the clusters that were condensed and a fasta file with the representative sequences of each cluster (but keeping all original genes, so no filtering is done).

# ## Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
from collections import Counter
from collections import defaultdict
import os
from scipy.stats import linregress
import matplotlib.pyplot as plt
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import shutil


def merge_clarc_results(paths_list, out_path, panaroo_true, ppanggo_true):

    all_core_clusters = pd.DataFrame()
    all_acc_clusters = pd.DataFrame()

    for path in paths_list:

        core_path = path+'/clarc_output/clarc_results/core_cluster_summary.csv'
        acc_path = path+'/clarc_output/clarc_results/accessory_cluster_summary.csv'

        current_core_df = pd.read_csv(core_path)
        cog_columns = [col for col in current_core_df.columns if col.startswith('COG')]
        current_core_df = current_core_df[cog_columns]

        current_acc_df = pd.read_csv(acc_path)
        cog_columns = [col for col in current_acc_df.columns if col.startswith('COG')]
        current_acc_df = current_acc_df[cog_columns]


        # Join dataframes with core clusters
        all_core_clusters = pd.concat([all_core_clusters, current_core_df], join='outer', ignore_index=True)
        all_acc_clusters = pd.concat([all_acc_clusters, current_acc_df], join='outer', ignore_index=True)

    # Join dataframe to include all clusters found by CLARC (core and accessory)
    all_clusters = pd.concat([all_core_clusters, all_acc_clusters], join='outer', ignore_index=True)

    # Check which clusters have COGs that appear in any other cluster and assign unique cluster group IDs

    def has_common_elements(row1, row2):
        set1 = set(row1.dropna())
        set2 = set(row2.dropna())
        return bool(set1 & set2)

    adj_list = defaultdict(list)
    for i in range(len(all_clusters)):
        for j in range(i + 1, len(all_clusters)):
            if has_common_elements(all_clusters.iloc[i], all_clusters.iloc[j]):
                adj_list[i].append(j)
                adj_list[j].append(i)

    def dfs(node, group_id, visited, group_ids):
        stack = [node]
        while stack:
            current = stack.pop()
            if not visited[current]:
                visited[current] = True
                group_ids[current] = group_id
                for neighbor in adj_list[current]:
                    if not visited[neighbor]:
                        stack.append(neighbor)

    visited = [False] * len(all_clusters)
    group_ids = [-1] * len(all_clusters)
    group_id = 0

    for i in range(len(all_clusters)):
        if not visited[i]:
            dfs(i, group_id, visited, group_ids)
            group_id += 1

    # Assign group IDs to the dataframe
    all_clusters['cluster_group'] = group_ids

    # Calculate the size of each group
    group_sizes = all_clusters['cluster_group'].value_counts().to_dict()

    # Map group sizes back to the dataframe
    all_clusters['cluster_group_size'] = all_clusters['cluster_group'].map(group_sizes)

    # Get clusters with COGs that only appear in one population
    clusters_all = all_clusters[all_clusters['cluster_group_size']==1]

    # Get clusters with COGs that only appear in multiple populations
    multpop_clusters = all_clusters[all_clusters['cluster_group_size']>1]

    grouped_clusters = multpop_clusters.groupby('cluster_group')

    count_included = 0
    count_ignored = 0
    cog_columns_long = [col for col in all_clusters.columns if col.startswith('COG')]
    longest_cluster_size = len(cog_columns_long)

    for group_name, group_df in grouped_clusters:

        # Get group dataframe
        x = group_df.reset_index(drop=True)

        ## Obtain row with longest group
        # Count the non-NA values in the first n columns for each row (this number depends on the longest cluster overall)
        non_na_counts = x.iloc[:, :longest_cluster_size].apply(lambda row: row.count(), axis=1)
        # Find the index of the row with the maximum non-NA count
        max_non_na_index = non_na_counts.idxmax()
        # Get the row with the most non-NA values
        row_with_most_elements = x.loc[max_non_na_index]
        # Get length of longest cluster
        longest_cluster_length = non_na_counts.max()

        ## Check the unique number of elements across all individual clusters
        # Extract the first n columns
        first_five_cols = x.iloc[:, :longest_cluster_size]
        # Flatten the values into a single list and remove NaN values
        all_elements = first_five_cols.values.flatten()
        all_elements = [element for element in all_elements if pd.notna(element)]
        # Get the unique elements
        unique_elements = set(all_elements)
        # Count the unique elements
        unique_count = len(unique_elements)

        ## Find clusters to keep and append them to the previous list of clusters to keep (the ones that appeared in only one population)

        if longest_cluster_length == unique_count:

            count_included = count_included+1

            # Convert the new row into a DataFrame
            row_with_most_elements_df = pd.DataFrame([row_with_most_elements])

            # Concatenate the new row to the existing DataFrame
            clusters_all = pd.concat([clusters_all, row_with_most_elements_df], ignore_index=True)

        else:
            count_ignored = count_ignored+1

    # Get list of unique COGs in the accessory clusters to keep

    # Extract the first n columns
    first_cols = clusters_all.iloc[:, :longest_cluster_size]

    # Flatten the values into a single list and remove NaN values
    all_elements = first_cols.values.flatten()
    all_elements = [element for element in all_elements if pd.notna(element)]

    # Get the unique elements
    all_unique = list(set(all_elements))

    ## Create legend for the clusters to condense

    conn_cluster_cogs_legend = pd.DataFrame()

    for col in clusters_all.columns:
        if 'COG' in col:
            # Add '-' after each COG value
            conn_cluster_cogs_legend[col] = clusters_all[col] + '-'

    # Concatenate COG values with '-' in between
    conn_cluster_cogs_legend['new_cluster_name_samecog'] = conn_cluster_cogs_legend.fillna('').sum(axis=1)

    for col in clusters_all.columns:
        if 'COG' in col:
            conn_cluster_cogs_legend[col] = clusters_all[col]

    ### Condense the identified CLARC clusters on the original Roary output

    # Define function to generate the condensed row, given a list of the COGs present in the cluster and the new cluster name
    # The new cluster name is input as the new 'Gene' name, the No. isolate and sequences are summed to obtain the new Avg seqs per isolate of the cluster
    # For all other columns, the values of each COG in the cluster are just concatenated

    def create_condensed_row(df, genes, new_gene_name):
        rows_to_combine = df[df['Gene'].isin(genes)]

        if rows_to_combine.empty:
            return None

        # Perform aggregation
        condensed_row = {
            'Gene': new_gene_name,
            'Non-unique Gene name': ', '.join(rows_to_combine['Non-unique Gene name'].fillna('').astype(str)),
            'Annotation': ', '.join(rows_to_combine['Annotation'].fillna('').astype(str)),
            'No. isolates': rows_to_combine['No. isolates'].sum(),
            'No. sequences': rows_to_combine['No. sequences'].sum(),
            'Avg sequences per isolate': rows_to_combine['No. sequences'].sum() / rows_to_combine['No. isolates'].sum() if rows_to_combine['No. isolates'].sum() != 0 else 0,
            **{col: ', '.join(rows_to_combine[col].fillna('').astype(str)) for col in df.columns[7:]}
        }

        return condensed_row

    def create_condensed_row_ppanggo(df, genes, new_gene_name):
        rows_to_combine = df[df['Gene'].isin(genes)]

        if rows_to_combine.empty:
            return None

        # Perform aggregation
        condensed_row = {
            'Gene': new_gene_name,
            **{col: min(rows_to_combine[col].fillna(0).sum(), 1) for col in df.columns[1:]}
        }

        return condensed_row

    # Define function to loop through the legend of clusters and condense the rows of each cluster in the original pangenome output

    def process_clusters(df, legend):

        legend_data = legend.to_dict(orient='records')
        gene_set = set(df['Gene'])

        condensed_rows = []
        gene_to_condensed = {}

        for row in legend_data:
            values = list(row.values())
            cogs_to_cluster = [x for x in values[:-1] if pd.notna(x)]
            new_cluster_name = row['new_cluster_name_samecog']

            # Skip if there are no COGs in the cluster (shouldn't happen but who knows)
            if not cogs_to_cluster:
                continue

            condensed_row = create_condensed_row(df, cogs_to_cluster, new_cluster_name)
            if condensed_row:
                condensed_rows.append(condensed_row)

                for gene in cogs_to_cluster:
                    gene_to_condensed[gene] = new_cluster_name


        filtered_df = df[~df['Gene'].isin(gene_to_condensed.keys())]
        cluster_new_rows = pd.DataFrame(condensed_rows)
        final_df = pd.concat([filtered_df, cluster_new_rows], ignore_index=True)

        return final_df

    def process_clusters_ppanggo(df, legend):

        legend_data = legend.to_dict(orient='records')
        gene_set = set(df['Gene'])

        condensed_rows = []
        gene_to_condensed = {}

        for row in legend_data:
            values = list(row.values())
            cogs_to_cluster = [x for x in values[:-1] if pd.notna(x)]
            new_cluster_name = row['new_cluster_name_samecog']

            # Skip if there are no COGs in the cluster (shouldn't happen but who knows)
            if not cogs_to_cluster:
                continue

            condensed_row = create_condensed_row_ppanggo(df, cogs_to_cluster, new_cluster_name)
            if condensed_row:
                condensed_rows.append(condensed_row)

                for gene in cogs_to_cluster:
                    gene_to_condensed[gene] = new_cluster_name


        filtered_df = df[~df['Gene'].isin(gene_to_condensed.keys())]
        cluster_new_rows = pd.DataFrame(condensed_rows)
        final_df = pd.concat([filtered_df, cluster_new_rows], ignore_index=True)

        return final_df

    # Now we want to drop the COGs in the clusters we just condensed and add the presence absence info of the newly defined clusters. The presence absence matrix selected will the the one in the data folder of the first path, but as a reminder, they should be the same in ALL paths.

    # Import original presence absence matrix results and filter them to generate cleaner presence absence matrix (but do NOT filter for core or accessory)
    # The filtering will depend on whether the input is from Roary or Panaroo

    if panaroo_true == 0 and ppanggo_true == 0:

        pres_abs_path = paths_list[0]+'/data/gene_presence_absence.csv'

        # Import output pres abs
        roary_all_c = pd.read_csv(pres_abs_path, low_memory=False)
        roary_all_c['Gene'] = roary_all_c['Gene'].str.replace(r"[ ,\'\"]", '_', regex=True)


        # First round of filtering by fragment length >150bp OR general isolate frequency >10%
        #tenp = (roary_all_c.shape[1]-14)/10 # Counting number of isolates (first 14 columns are metadata)
        #roary_onefilt = roary_all_c[(roary_all_c['Avg group size nuc'] >= 150) | (roary_all_c['No. isolates'] >= tenp)]
        roary_onefilt = roary_all_c.copy()

        # Now make gene names the indeces
        roary_onefilt.set_index('Gene', inplace=True)

        # Drop all columns that are not an isolate
        roary_isol = roary_onefilt.iloc[:,13:]

        # Now replace NaN values for 0 and any other value for 1
        roary_isol[~roary_isol.isnull()] = 1
        roary_isol[roary_isol.isnull()] = 0

        # Switch rows to columns for ease
        pres_abs_isol = roary_isol.transpose()
        pres_abs_isol.index.name='Accession'

        # Condense clusters on original Roary output, using pre-defined functions
        roary_onefilt_condensed = roary_onefilt.copy().reset_index(drop=False)
        roary_onefilt_condensed = process_clusters(roary_onefilt_condensed, conn_cluster_cogs_legend)

        clarc_merge_path = out_path+"/clarc_merge_results"

        # Check if the directory already exists (create if it doesn't)
        if not os.path.exists(clarc_merge_path):
            os.makedirs(clarc_merge_path)

        # Export results
        roary_og_merge_clarced_path = clarc_merge_path+'/gene_presence_absence_clarc_merged.csv'
        roary_onefilt_condensed.to_csv(roary_og_merge_clarced_path, index=False)

    elif panaroo_true == 1:

        pres_abs_path = paths_list[0]+'/data/gene_presence_absence_roary.csv'

        # Import and filter panaroo results
        panaroo_all_c = pd.read_csv(pres_abs_path, low_memory=False)
        panaroo_all_c['Gene'] = panaroo_all_c['Gene'].str.replace(r"[ ,\'\"]", '_', regex=True)

        # First round of filtering by fragment length >150bp OR general isolate frequency >10%
        #tenp = (panaroo_all_c.shape[1]-14)/10 # Counting number of isolates (first 14 columns are metadata)
        #panaroo_onefilt = panaroo_all_c[(panaroo_all_c['Avg group size nuc'] >= 150) | (panaroo_all_c['No. isolates'] >= tenp)]
        panaroo_onefilt = panaroo_all_c.copy()

        # Now make gene names the indeces
        panaroo_onefilt.set_index('Gene', inplace=True)

        # Drop all columns that are not an isolate
        panaroo_isol = panaroo_onefilt.iloc[:,13:]

        # Now replace NaN values for 0 and any other value for 1
        panaroo_isol[~panaroo_isol.isnull()] = 1
        panaroo_isol[panaroo_isol.isnull()] = 0

        # Switch rows to columns for ease
        pres_abs_isol = panaroo_isol.transpose()
        pres_abs_isol.index.name='Accession'

        # Condense clusters on original Roary output, using pre-defined functions
        panaroo_onefilt_condensed = panaroo_onefilt.copy().reset_index(drop=False)
        panaroo_onefilt_condensed = process_clusters(panaroo_onefilt_condensed, conn_cluster_cogs_legend)

        clarc_merge_path = out_path+"/clarc_merge_results"

        # Check if the directory already exists (create if it doesn't)
        if not os.path.exists(clarc_merge_path):
            os.makedirs(clarc_merge_path)

        # Export results
        panaroo_og_merge_clarced_path = clarc_merge_path+'/gene_presence_absence_clarc_merged.csv'
        panaroo_onefilt_condensed.to_csv(panaroo_og_merge_clarced_path, index=False)

    elif ppanggo_true == 1:

        # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)
        pres_abs_path = paths_list[0]+'/ppanggo_data/gene_presence_absence.Rtab'

        # Import Ppanggo output
        igopan_all_ppanggo = pd.read_csv(pres_abs_path, sep='\t', comment='#')

        # Bakta often has spaces and special characters in its gene names, which can cause problems in the downstream analyses. So here we can turn them into underscores. I'm still deciding if this update is necessary. I think as long as the fasta files pick up the whole name, it should be fine.
        # Update: this was absolutely necessary
        igopan_all_ppanggo["Gene"] = igopan_all_ppanggo["Gene"].str.replace("[ ,\'\"]", "_", regex=True)

        # Get list of PPanGGolin output names in a list
        panppanggo_ids_list =  list(igopan_all_ppanggo["Gene"])

        # Now make gene names the indeces
        igopan_all_ppanggo.set_index('Gene', inplace=True)

        # Switch rows to columns for ease
        pres_abs_isol = igopan_all_ppanggo.transpose()
        pres_abs_isol.index.name='Accession'

        # Condense clusters on PPanGGolin output, using pre-defined functions
        ppanggo_onefilt_condensed = igopan_all_ppanggo.copy().reset_index(drop=False)
        ppanggo_onefilt_condensed = process_clusters_ppanggo(ppanggo_onefilt_condensed, conn_cluster_cogs_legend)

        clarc_merge_path = out_path+"/clarc_merge_results"

        # Check if the directory already exists (create if it doesn't)
        if not os.path.exists(clarc_merge_path):
            os.makedirs(clarc_merge_path)

        # Export results
        ppanggo_og_merge_clarced_path = clarc_merge_path+'/gene_presence_absence_clarc_merged.csv'
        ppanggo_onefilt_condensed.to_csv(ppanggo_og_merge_clarced_path, index=False)

    # Add empty column to dataframe with cog pairs that exclude each other
    pres_abs_newclusters = pd.DataFrame()

    for index, row in conn_cluster_cogs_legend.iterrows():

        # Get COG names for row and new name
        cog_columns = [f'COG{i}' for i in range(1, longest_cluster_size+1)]
        cogs = [row[cog] for cog in cog_columns if pd.notnull(row[cog])]
        new_name = row['new_cluster_name_samecog']

     # Initialize dataframe
        df = pres_abs_isol[cogs].copy()  # Use .copy() to avoid SettingWithCopyWarning

        # Concatenate COGs
        df[new_name] = df.sum(axis=1)
        pres_abs_newclusters = pd.concat([pres_abs_newclusters, df[new_name]], axis=1)
        # Since only the linkage of a specific subpopulation might be being used, if the the COGs coexist in any other sample then they are counted as present
        pres_abs_newclusters[pres_abs_newclusters > 1] = 1

    # Drop cluster COGs from presence absence matrix
    nonredundant_presabs = pres_abs_isol.drop(columns=all_unique)

    # Now append new groups
    samecog_clustered_presabs = pd.concat([nonredundant_presabs, pres_abs_newclusters], axis=1)

    ## Save the summary output files for the merged clusters
    out_pres_abs = clarc_merge_path+'/presence_absence_clarc_merged_binary.csv'
    samecog_clustered_presabs.to_csv(out_pres_abs, index=True)

    out_text = clarc_merge_path+'/clarc_merge_cluster_summary.txt'
    with open(out_text, "a") as file:
        file.write(f'Total COG clusters identified and merged by CLARC: '+f'{len(clusters_all)}\n')
        file.write(f'Unique COGs in the clusters: '+f'{len(all_unique)}\n')

    all_cluster_path = out_path+"/clarc_merge_results/merge_cluster_summary.csv"
    clusters_all.to_csv(all_cluster_path, index=False)

    cog_name_list = list(all_unique)

    # Define the file path
    list_path = out_path+"/clarc_merge_results/list_cluster_cogs.txt"

    # Open the file in write mode
    with open(list_path, "w") as file:
        # Write each element of the list to the file
        for item in cog_name_list:
            file.write(item + "\n")  # Add a newline after each item

    ## Get condensed fasta file

    # Get fasta file with the representative sequences of the new accessory gene definitions
    # Get dataframe with all accessory COGs and their length

    if ppanggo_true == 0:

        original_fasta_path = paths_list[0]+'/data/pan_genome_reference.fa'

    elif ppanggo_true == 1:

        original_fasta_path = paths_list[0]+'/ppanggo_data/all_nucleotide_families.fasta'


    def fasta_coglen_dataframe(fasta_file):

        seq_ids = []
        seq_lengths = []

        for record in SeqIO.parse(fasta_file, "fasta"):
            entry_description = record.description
            parts = entry_description.split(' ', 1)

            # If there is no description after the ID, just use the ID
            if len(parts) == 1:
                cog_name = parts[0]
            else:
                cog_name = parts[1]

            cog_name_fixed = (cog_name.replace(' ', '_').replace(',', '_').replace("'", '_').replace('"', '_'))

            seq_ids.append(cog_name_fixed)
            seq_lengths.append(len(record.seq))

        df = pd.DataFrame({"cog_name": seq_ids, "length": seq_lengths})

        return df

    # Get dataframe with length of each COG evaluated by CLARC (previous accessory COGs)

    # Parse the fasta file and create the DataFrame
    cog_len = fasta_coglen_dataframe(original_fasta_path)

    conn_cluster_cogs_legend = conn_cluster_cogs_legend.reset_index(drop=True)

    # Identify longest COG in each CLARC cluster

    for index, row in conn_cluster_cogs_legend.iterrows():
        longest_cog = ''
        max_len = 0

        for i in range(1, longest_cluster_size + 1):
            cog = row[f'COG{i}']

            if not pd.isna(cog):  # Check if cog is not NaN
                len_row = cog_len.query(f"cog_name == '{cog}'")
                if not len_row.empty:
                    leng = len_row['length'].values[0]

                    if leng > max_len:
                        longest_cog = cog
                        max_len = leng

        conn_cluster_cogs_legend.at[index, 'longest_cog'] = longest_cog

    # Get important lists
    longest_cogs = list(conn_cluster_cogs_legend['longest_cog'])

    # Get fasta file with the representative sequences for the CLARC COG definitions - accessory and core

    ## For accessory:
    # The COGs identified as part of core clusters will be dropped
    # The COGs identified as part of an accessory cluster, but not the longest cog will be dropped
    # The COGs identified as part of an accessory cluster and representing the longest COG will have their name reannotated with the cluster name

    # Important: Bakta allows spaces and commas in their COG names, which messes with the way in which biopython reads the record id
    # So, all spaces and commas in COG names are replaced by underscores throughout the CLARC pipeline

    clarc_fasta_path_acc = out_path+'/clarc_merge_results/pan_genome_reference_clarc_merge.fasta'

    fin = open(original_fasta_path, 'r')
    fout = open(clarc_fasta_path_acc, 'w')

    for record in SeqIO.parse(fin,'fasta'):

        entry_description = record.description
        parts = entry_description.split(' ', 1)

        if len(parts) == 1:
            cog_name = parts[0]
        else:
            cog_name = parts[1]

        cog_name_fixed = cog_name.replace(' ', '_').replace(',', '_').replace("'", '_').replace('"', '_')

        if cog_name_fixed in all_unique:

            if cog_name_fixed in longest_cogs:

                x = conn_cluster_cogs_legend.query(f" longest_cog == '{cog_name_fixed}'")
                cluster_name = x['new_cluster_name_samecog'].values[0]
                record.id = cluster_name
                SeqIO.write(record, fout, 'fasta')

        else:
            record.id = cog_name_fixed
            #record.description = ' '
            SeqIO.write(record, fout, 'fasta')

    fin.close()
    fout.close()

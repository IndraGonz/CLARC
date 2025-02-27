#!/usr/bin/env python
# coding: utf-8

## Condensing COG definitions with CLARC method

# Here we use linkage, functional annotations and sequence identity to identify clusters of COGs that are the same gene.

# ## Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
import itertools
import fastcluster
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import shutil
import math
from collections import Counter, defaultdict
import time

# Starting from CLARC 1.1.3
# Revision list:
#  - We don't use or need the log co-occur matrix
#- Made various parts of the code more efficient including:
  #- In the step where the blastn results are being integrated into the pairwise P11 dataframe, I've optimized the way the query is done. It was previously done very inefficiently with itterrows, but now I build a lookup table which really improves performance. This was a big bottleneck.
#- Added a BLASTN e-value minimum for the connections (<1e-10)
#- Added a user defined threshold where the linkage parameter can be varied. For example, instead of finding pairs that perfectly exclude each other, you can find pairs that co-occur <1% of the time.
#- Added a user defined threshold to relax the graph connectivity threshold. Instead of selecting fully connected clusters, now the user can select clusters that are 90% connected for example. However, this might make it so that some COGs appear in multiple clusters. To address this, we identify any COGs that appear in multiple clusters, and only keep them in the most connected of the clusters.
#- Added option to process PPanGGolin output

# ## Identify CLARC clusters, condense them and generate output files

def clarc_cleaning(in_path, out_path, panaroo_true, ppanggo_true, acc_upper, acc_lower, core_lower, clarc_identity, linkage_cut, connection_cut):

    # Import counts for COG pairs that co-occur
    linkage_p11_path = out_path+"/linkage/acc_p11_matrix.csv"
    linkage_p11_result = pd.read_csv(linkage_p11_path ,low_memory=False, index_col=0)

    # Check COG pairs that never co-occur
    p11_coglist = linkage_p11_result.melt(ignore_index = False)
    p11_coglist = p11_coglist.reset_index()
    p11_coglist.columns = ['COG1','COG2','p11']

    # Import presence absence matrix
    pres_abs_path = out_path+'/population_accessory_presence_absence.csv'
    pan_acc = pd.read_csv(pres_abs_path, index_col=0, low_memory=False)

    # Get number of co-occurrence events that would represent the percentage threshold given the number of samples
    num_samples = len(pan_acc)
    coccur_events = math.ceil(num_samples * linkage_cut)

    # Get list of pairs that co-occur in less than the specified linkage threshold
    p11_coglist_zeros = p11_coglist[p11_coglist["p11"] <= coccur_events]

    # Now we have to drop the duplicates from these dataframes

    # p11
    df_p11 = pd.DataFrame(np.sort(p11_coglist_zeros[['COG1','COG2']].values,1), columns=['COG1', 'COG2']).drop_duplicates()
    p11_coglist_zeros_unique = p11_coglist_zeros.merge(df_p11, on=['COG1', 'COG2'])

    # Calculate COG frequencies in the population

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
    blast_path = out_path + '/acc_blastn/blastn_acccogs_allvall.tsv'
    blastn_acccogs = pd.read_table(blast_path, header=None)

    # TSV output doesn't have column labels, so here I add them
    column_labels = ['query_seq_ID', 'subject_seq_ID', 'percentage_identical_matches', 'align_length',
                     'num_mismatches', 'gap_open', 'align_start_query', 'align_end_query',
                     'align_start_subject', 'align_end_subject', 'e-value', 'bit_score']
    blastn_acccogs.columns = column_labels

    # Resolve duplicates: Keep the row with the **largest align_length** for each (query_seq_ID, subject_seq_ID) pair
    blastn_acccogs = blastn_acccogs.sort_values('align_length', ascending=False).drop_duplicates(subset=['query_seq_ID', 'subject_seq_ID'])

    # Create lookup tables using MultiIndex
    blast_lookup = blastn_acccogs.set_index(['query_seq_ID', 'subject_seq_ID'])[['percentage_identical_matches', 'align_length', 'e-value']]

    # Convert tuples (COG1, COG2) into an index for lookup
    lookup_keys = list(zip(p11_coglist_freq['COG1'], p11_coglist_freq['COG2']))

    # Add new columns to p11_coglist_freq
    p11_coglist_freq[['blastn_iden_per', 'blastn_align_length', 'blastn_e_value']] = blast_lookup.reindex(lookup_keys).values

    # Get length of each queried acc COG

    # Extract unique COG lengths where COG1 == COG2
    cog_lengths = p11_coglist_freq[p11_coglist_freq['funct_group_cog1'] == p11_coglist_freq['funct_group_cog2']]
    cog_length_dict = dict(zip(cog_lengths['funct_group_cog1'], cog_lengths['blastn_align_length']))

    # Map lengths to COG1 and COG2 to create new columns
    p11_coglist_freq['cog_len_1'] = p11_coglist_freq['funct_group_cog1'].map(cog_length_dict)
    p11_coglist_freq['cog_len_2'] = p11_coglist_freq['funct_group_cog2'].map(cog_length_dict)

    # Same eggnog functional group
    samecog_condense = p11_coglist_freq[p11_coglist_freq['eggnog_same_flag'] == 1]

    # High blastn sequence identity
    samecog_condense = samecog_condense[samecog_condense['blastn_iden_per'] >= clarc_identity]

    # Drop rows where COG1 is equal to COG2
    samecog_condense = samecog_condense[samecog_condense['COG1'] != samecog_condense['COG2']]

    # Keep only the entries where the alignment length covers at least 80% of the length of the shortest COG
    #samecog_condense = samecog_condense[(samecog_condense["blastn_align_length"] / samecog_condense[["cog_len_1", "cog_len_2"]].min(axis=1)) > 0.80]

    # Keep only entries where the e-value is above 1e-10
    samecog_condense = samecog_condense[samecog_condense['blastn_e_value'] <= 1e-10]

    # Find connected clusters above a connectivity threshold specified by the user (default is 100%, meaning fully connected)
    p11_zero_cog1_list = list(samecog_condense['COG1'])
    p11_zero_cog2_list = list(samecog_condense['COG2'])

    samegene_cogs = p11_zero_cog1_list + p11_zero_cog2_list
    samegene_cogs_list = list(set(samegene_cogs))  # Unique COGs

    # Initialize DataFrame to store connected groups
    samegene_connected_cogs = pd.DataFrame(columns=['connected_COGs', 'group_number', 'group_length'])
    group = 0

    # Iterate through each unique COG
    for cog in samegene_cogs_list:
        # Get all COGs connected to `cog`
        con = samecog_condense.query(f"COG1 == '{cog}' or COG2 == '{cog}'")
        x = sorted(pd.concat([con['COG1'], con['COG2']]).unique())
        xl = x.copy()
        xl.remove(cog)  # Remove itself from the list

        # Ensure full connectivity or threshold-based grouping
        flag = 0

        for i in xl:
            coni = samecog_condense.query(f"COG1 == '{i}' or COG2 == '{i}'")
            y = sorted(pd.concat([coni['COG1'], coni['COG2']]).unique())

            # If threshold = 1.0, enforce strict full connectivity (original logic)
            if connection_cut == 1.0:
                if set(y) == set(x):  # Must match exactly (fully connected)
                    flag += 1
            else:
                # If threshold < 1.0, allow slightly looser grouping
                common_connections = len(set(y) & set(x))  # Number of shared connections
                total_possible_connections = len(x)  # Expected connections
                connection_ratio = common_connections / total_possible_connections
                if connection_ratio >= connection_cut:
                    flag += 1

        # Only add as a group if all members meet the threshold requirement
        if flag == len(xl) if connection_cut == 1.0 else flag / len(xl) >= connection_cut:
            group += 1
            samegene_connected_cogs.loc[group-1, 'connected_COGs'] = x
            samegene_connected_cogs.loc[group-1, 'group_number'] = group
            samegene_connected_cogs.loc[group-1, 'group_length'] = len(x)

    # Convert the connected_COGs column to tuples for deduplication
    samegene_connected_cogs['connected_COGs'] = samegene_connected_cogs['connected_COGs'].apply(tuple)

    # Remove duplicate clusters
    samegene_connected_cogs = samegene_connected_cogs.drop_duplicates(subset='connected_COGs', keep="first")

    # Sort by largest group first
    samegene_connected_cogs = samegene_connected_cogs.sort_values(by=['group_length'], ascending=False)

    # Reset index
    samegene_connected_cogs = samegene_connected_cogs.reset_index(drop=True)

    ### STEP 1: FIND COGs APPEARING IN MULTIPLE GROUPS ###
    # Since we are relaxing the connectivity parameter, this might happen
    # Flatten the list of all COGs
    all_cogs = [cog for group in samegene_connected_cogs["connected_COGs"] for cog in group]

    # Count occurrences of each COG
    cog_counts = Counter(all_cogs)

    # Find COGs that appear in multiple groups
    duplicate_cogs = {cog: count for cog, count in cog_counts.items() if count > 1}

    ### STEP 2: FIX DUPLICATE COGs (Assign to Most Connected Group) ###
    # Create a dictionary to track which COG belongs to which group based on **highest number of direct connections**
    cog_to_group = {}

    # Track COG connections within each group
    cog_connection_counts = defaultdict(lambda: defaultdict(int))

    # Count how many connections each COG has in each group
    for index, row in samegene_connected_cogs.iterrows():
        group_number = row['group_number']
        connected_cogs = row['connected_COGs']

        for cog in connected_cogs:
            for other_cog in connected_cogs:
                if cog != other_cog:
                    cog_connection_counts[cog][group_number] += 1  # Count number of direct links in each group

    # Assign each COG to the **group where it has the most direct connections**
    for cog, group_dict in cog_connection_counts.items():
        best_group = max(group_dict, key=group_dict.get)  # Group with highest connections
        cog_to_group[cog] = best_group

    # Rebuild clusters ensuring each COG appears only once, keeping them in the **most connected group**
    unique_clusters = defaultdict(list)

    for cog, group in cog_to_group.items():
        unique_clusters[group].append(cog)

    # Convert back to DataFrame
    new_cluster_list = [{'connected_COGs': sorted(v), 'group_length': len(v)} for k, v in unique_clusters.items()]
    samegene_connected_cogs = pd.DataFrame(new_cluster_list)

    if samegene_connected_cogs.empty:

            clarc_summary_path = out_path+"/clarc_results"

            if not os.path.exists(clarc_summary_path):
                os.makedirs(clarc_summary_path)

            out_text = clarc_summary_path+'/clarc_cluster_summary.txt'
            with open(out_text, "a") as file:
                file.write(f'Total COG clusters identified by CLARC: '+f'{0}\n')
                file.write(f'Core COG clusters: '+f'{0}\n')
                file.write(f'Unique COGs in core clusters: '+f'{0}\n')
                file.write(f'Accessory COG clusters: '+f'{0}\n')
                file.write(f'Unique COGs in accessory clusters: '+f'{0}\n')

            sys.exit("No CLARC clusters found")

    else:

        # Sort by largest group first
        samegene_connected_cogs = samegene_connected_cogs.sort_values(by=['group_length'], ascending=False).reset_index(drop=True)

        # Remove clusters with only one COG
        samegene_connected_cogs = samegene_connected_cogs[samegene_connected_cogs['group_length'] > 1]

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

        cog_columns = [f"COG{i}" for i in range(1, max_clusters + 1)]

        cluster_cogs_list = []
        for col in cog_columns:
            cog_list_a = list(samegene_connected_cogs_exp[col])
            cog_list = list(filter(None, cog_list_a))
            cluster_cogs_list.extend(cog_list)

        cluster_cogs_unique = set(cluster_cogs_list)

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
            file.write(f'Total COG clusters identified by CLARC: '+f'{len(samegene_connected_cogs_exp)}\n')
            file.write(f'Core COG clusters: '+f'{len(conn_cluster_cogs_core)}\n')
            file.write(f'Unique COGs in core clusters: '+f'{len(cluster_corecogs_unique)}\n')
            file.write(f'Accessory COG clusters: '+f'{len(conn_cluster_cogs_acc)}\n')
            file.write(f'Unique COGs in accessory clusters: '+f'{len(cluster_acccogs_unique)}\n')

        core_cluster_path = out_path+"/clarc_results/core_cluster_summary.csv"
        conn_cluster_cogs_core.to_csv(core_cluster_path, index=False)

        acc_cluster_path = out_path+"/clarc_results/accessory_cluster_summary.csv"
        conn_cluster_cogs_acc.to_csv(acc_cluster_path, index=False)

        # Condense all identified clusters in original pangenome results - New presence absence matrix

        # Get legend of new clusters to condense

        conn_cluster_cogs_legend = pd.DataFrame()
        for col in samegene_connected_cogs_exp.columns:
            if 'freq' not in col:
                # Add '-' after each COG value
                conn_cluster_cogs_legend[col] = samegene_connected_cogs_exp[col] + '-'

       # Concatenate COG values with '-' in between
        conn_cluster_cogs_legend['new_cluster_name_samecog'] = conn_cluster_cogs_legend.fillna('').sum(axis=1)

        for col in samegene_connected_cogs_exp.columns:
            if 'freq' not in col:
                conn_cluster_cogs_legend[col] = samegene_connected_cogs_exp[col]

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

        # Now we want to drop the COGs in the clusters we just condensed and add the presence absence info of the newly defined clusters

        # Import original presence absence matrix results and filter them to generate cleaner presence absence matrix (but do NOT filter for core or accessory)
        # The filtering will depend on whether the input is from Roary or Panaroo

        if panaroo_true == 0 and ppanggo_true == 0:

            pres_abs_path = in_path+'/gene_presence_absence.csv'

            # Import output pres abs
            roary_all_c = pd.read_csv(pres_abs_path, low_memory=False)
            roary_all_c['Gene'] = roary_all_c['Gene'].str.replace(r"[ ,\'\"]", '_', regex=True)

            # Get list of Roary output names in a list
            panroary_ids_list =  list(roary_all_c["Gene"])

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

            # Export results
            roary_og_clarced_path = out_path+'/clarc_results/roary_clarc_gene_presence_absence.csv'
            roary_onefilt_condensed.to_csv(roary_og_clarced_path, index=False)

        elif panaroo_true == 1:

            pres_abs_path = in_path+'/gene_presence_absence_roary.csv'

            # Import and filter panaroo results
            panaroo_all_c = pd.read_csv(pres_abs_path, low_memory=False)
            panaroo_all_c['Gene'] = panaroo_all_c['Gene'].str.replace(r"[ ,\'\"]", '_', regex=True)


            # Get list of Roary output names in a list
            panaroo_ids_list =  list(panaroo_all_c["Gene"])

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

            # Export results
            panaroo_og_clarced_path = out_path+'/clarc_results/panaroo_clarc_gene_presence_absence.csv'
            panaroo_onefilt_condensed.to_csv(panaroo_og_clarced_path, index=False)

        elif ppanggo_true == 1:

            # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)
            pres_abs_path = in_path+"/gene_presence_absence.Rtab"

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

            # Export results
            ppanggo_og_clarced_path = out_path+'/clarc_results/ppanggolin_clarc_gene_presence_absence.csv'
            ppanggo_onefilt_condensed.to_csv(ppanggo_og_clarced_path, index=False)

        # Add empty column to dataframe with cog pairs that exclude each other
        pres_abs_newclusters = pd.DataFrame()

        for index, row in conn_cluster_cogs_legend.iterrows():

            # Get COG names for row and new name
            cog_columns = [f'COG{i}' for i in range(1, max_clusters+1)]
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
        nonredundant_presabs = pres_abs_isol.drop(columns=cluster_cogs_unique)

        # Now append new groups
        samecog_clustered_presabs = pd.concat([nonredundant_presabs, pres_abs_newclusters], axis=1)

        # We save the output
        clarc_all_condense_path = out_path+'/clarc_results/clarc_condensed_gene_pres_abs_binary.csv'
        samecog_clustered_presabs.to_csv(clarc_all_condense_path, index=True)

        # Awesome, now we can filter again for accessory and core genes

        sample_needed_path = in_path+'/needed_sample_names.txt'

        # Get list of sample names in the subpopulation to analyse
        with open(sample_needed_path, 'r') as file:
            acc_needed_list = file.read().splitlines()

        # Get only isolates in the list
        genefreq_mat_filt = samecog_clustered_presabs[samecog_clustered_presabs.index.isin(acc_needed_list)]

        # Now we filter by frequency per each dataset (this can change depending on the datasets we want to include)

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

        freq_cog = get_freq(genefreq_mat_filt)

        # Now filter to only keep cogs with frequency between set thresholds, default between 5 and 95% in the given sample subset
        acc_cog = freq_cog.loc[(freq_cog['freq'] >= acc_lower) & (freq_cog['freq'] <= acc_upper)]

        # Now filter to only keep cogs with frequency over set threshold, default over 95% in the given sample subset
        core_cog = freq_cog.loc[(freq_cog['freq'] > core_lower)]

        # Get the number of accessory cogs per dataset
        acccog_num = acc_cog.shape[0]
        corecog_num = core_cog.shape[0]

        # Get accessory cog names as lists
        acccog_name_list = list(acc_cog["COG_name"])
        corecog_name_list = list(core_cog["COG_name"])

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

        # Only get presence absence matrix for accessory genes
        genefreq_meta_filt_acc = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(acccog_name_list)]]
        genefreq_meta_filt_core = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(corecog_name_list)]]

        ## Export results to output folder
        acc_path = out_path+'/clarc_results/clarc_population_acc_presence_absence.csv'
        core_path = out_path+'/clarc_results/clarc_population_core_presence_absence.csv'

        genefreq_meta_filt_acc.to_csv(acc_path, index=True)
        genefreq_meta_filt_core.to_csv(core_path, index=True)

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

        # Get dataframe with length of each COG evaluated by CLARC (previous accessory COGs)

        # Provide the path to your fasta file
        original_fasta_path = out_path+'/accessory_rep_seqs.fasta'

        # Parse the fasta file and create the DataFrame
        cog_len = fasta_coglen_dataframe(original_fasta_path)

        conn_cluster_cogs_legend = conn_cluster_cogs_legend.reset_index(drop=True)

        # Identify longest COG in each CLARC cluster

        for index, row in conn_cluster_cogs_legend.iterrows():

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

                conn_cluster_cogs_legend.at[index, 'longest_cog'] = longest_cog

        # Get important lists
        longest_cogs = list(conn_cluster_cogs_legend['longest_cog'])

        # Get fasta file with the representative sequences for the CLARC COG definitions - accessory and core

        ## For accessory:
        # The COGs identified as part of core clusters will be dropped
        # The COGs identified as part of an accessory cluster, but not the longest cog will be dropped
        # The COGs identified as part of an accessory cluster and representing the longest COG will have their name reannotated with the cluster name

        clarc_fasta_path_acc = out_path+'/clarc_results/clarc_acc_cog_seqs.fasta'

        fin = open(original_fasta_path, 'r')
        fout = open(clarc_fasta_path_acc, 'w')

        for record in SeqIO.parse(fin,'fasta'):

            cog_name = record.id

            if cog_name not in cluster_corecogs_unique:

                if cog_name in cluster_acccogs_unique:

                    if cog_name in longest_cogs:

                        x = conn_cluster_cogs_legend.query(f" longest_cog == '{cog_name}'")
                        cluster_name = x['new_cluster_name_samecog'].values[0]
                        record.id = cluster_name
                        SeqIO.write(record, fout, 'fasta')

                else:
                    SeqIO.write(record, fout, 'fasta')

        fin.close()
        fout.close()

        ## For core:
        # The COGs identified as part of an  cluster, but not the longest cog will be dropped
        # The COGs identified as part of a core cluster and representing the longest COG in their cluster will be appended to the fasta file with the core sequences, with its names updated to represent the cluster name

        clarc_fasta_path_core = out_path+'/clarc_results/clarc_core_cog_seqs.fasta'
        original_fasta_path_core = out_path+'/core_rep_seqs.fasta'

        # Copy current core sequences into new output file
        shutil.copyfile(original_fasta_path_core, clarc_fasta_path_core)

        fin = open(original_fasta_path, 'r')
        fout = open(clarc_fasta_path_core, 'a') # Append mode, according to google

        for record in SeqIO.parse(fin,'fasta'):

            cog_name = record.id

            if cog_name in cluster_corecogs_unique:

                if cog_name in longest_cogs:

                    x = conn_cluster_cogs_legend.query(f" longest_cog == '{cog_name}'")
                    cluster_name = x['new_cluster_name_samecog'].values[0]
                    record.id = cluster_name
                    SeqIO.write(record, fout, 'fasta')

        fin.close()
        fout.close()

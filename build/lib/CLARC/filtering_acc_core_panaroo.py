#!/usr/bin/env python
# coding: utf-8

# # Filtering pangenome presence absence matrix to get accessory and core genes from Panaroo output

# This code has a function that takes the 'gene_presence_absence.csv' output from Panaroo and filters to obtain a presence absence matrix of only the accessory genes of a population of samples, where the columns are the accessory genes and the indeces are the accession numbers of the samples selected for the analysis. It reads a text file containing the sample names that are to be selected. By default, accessory genes are defined as those in 5-95% of the selected samples. However, these threshholds can be changed by specifying when running CLARC.
#
#
# Additionally, from the 'pan_genome_reference.fa' file it creates a fasta file with the representative sequences of only the accessory genes. But...
#
# Unlike Roary, Panaroo does not keep a representative for _all_ genes in its 'pan_genome_reference.fa' file. So, this filtering script requires an additional input than the Roary one which is the 'gene_data.csv' Panaroo output. Using these inputs, this filtering script finds the length of the longest sequence of the COGs that were not found in the 'pan_genome_reference.fa' file and the sample it came from. From there is uses the 'gene_data.csv' file to obtain the representative sequence. This additional process is computationally intensive due to the size of the gene_data file and the iterative searched. Thus, this process is parallelized. Still, the time to run the filtering script for Panaroo will probably take significantly longer than for Roary.
#
# There is also a separate function that does all of this but for the core genes, which by default are defined as those in >95% of the samples in the population.

# ### Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import os

# ### Define functions

def get_pop_acc_pres_abs_panaroo(data_path, out_path, acc_upper, acc_lower):

    pres_abs_path = data_path+"/gene_presence_absence_roary.csv"
    sample_needed_path = data_path+'/needed_sample_names.txt'

    # Import Panaroo output
    igopan_all_panaroo = pd.read_csv(pres_abs_path, low_memory=False)

    # Get list of panaroo output names in a list
    panpanaroo_ids_list =  list(igopan_all_panaroo["Gene"])

    # First round of filtering by fragment length >150bp OR general isolate frequency >10%
    tenp = (igopan_all_panaroo.shape[1]-14)/10 # Counting number of isolates (first 14 columns are metadata)

    panaroo_onefilt = igopan_all_panaroo[(igopan_all_panaroo['Avg group size nuc'] >= 150) | (igopan_all_panaroo['No. isolates'] >= tenp)]

    # Now make gene names the indeces
    panaroo_onefilt.set_index('Gene', inplace=True)

    # Drop all columns that are not an isolate
    panaroo_isol = panaroo_onefilt.iloc[:,13:]

    # Now replace NaN values for 0 and any other value for 1
    panaroo_isol[~panaroo_isol.isnull()] = 1
    panaroo_isol[panaroo_isol.isnull()] = 0

    # Switch rows to columns for ease
    panaroo_genefreq_matrix = panaroo_isol.transpose()
    #panaroo_genefreq_matrix.columns = panaroo_genefreq_matrix.columns.str.replace('~~~', '-')
    panaroo_genefreq_matrix.index.name='Accession'

    # #Get list of sample names in the subpopulation to analyse
    with open(sample_needed_path, 'r') as file:
        acc_needed_list = file.read().splitlines()

    # Get only isolates in the list
    genefreq_mat_filt = panaroo_genefreq_matrix[panaroo_genefreq_matrix.index.isin(acc_needed_list)]

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

    # population all
    freq_cog = get_freq(genefreq_mat_filt)

    # Now filter to only keep cogs with frequency in a given range within the given sample subset
    acc_cog = freq_cog.loc[(freq_cog['freq'] >= acc_lower) & (freq_cog['freq'] <= acc_upper)]

    # Get the number of accessory cogs per dataset
    acccog_num = acc_cog.shape[0]

    # Get accessory cog names as lists
    acccog_name_list = list(acc_cog["COG_name"])

    # Only get presence absence matrix for accessory genes
    genefreq_meta_filt_acc = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(acccog_name_list)]]

    ## Export results to output folder
    csv_path = out_path+f'/population_accessory_presence_absence.csv'

    genefreq_meta_filt_acc.to_csv(csv_path, index=True)

    ### Now get fasta file with only accessory genes

    # Define the input and output file paths
    seq_path = data_path+"/pan_genome_reference.fa"
    acc_fasta_path = out_path+f"/accessory_rep_seqs.fasta"

    cog_id_list = []

    # Open the input and output files
    with open(seq_path, "r") as infile, open(acc_fasta_path, "w") as outfile:
        # Loop through each record in the input FASTA file
        for record in SeqIO.parse(infile, "fasta"):
            # Check if the sequence ID is in the list of sequence IDs to keep
            cog_id = record.description
            if cog_id in acccog_name_list:
                # Write the record to the output file
                record.id = cog_id
                record.description = ''
                SeqIO.write(record, outfile, "fasta")
                cog_id_list.append(cog_id)

    ### Find the rep sequences of the genes not represented in the pan_genome_reference.fa Panaroo output

    # Function to write fasta files
    def fasta_out(rec):
        return(">"+str(rec.description)+"\n"+str(rec.seq))

    # Function to format the header for comparison
    def format(rec):
        return(rec.id)

    # Function to extract sequence identifiers from a fasta file
    def extract_sequence_identifiers(fasta_file_path):
        identifiers = []
        with open(fasta_file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    # Extract the sequence identifier (remove the ">")
                    identifier = line[1:]
                    identifiers.append(identifier)
        return identifiers

    # Function to find the representative sequence
    def find_rep_seq(file_path, queries, max_len_seq):
        chunk_iter = pd.read_csv(file_path, chunksize=1000)
        max_sequence_length = 0
        max_sequence_row = None

        # Loop through each chunk
        for chunk in chunk_iter:
            # Filter rows based on 'gff_file' and 'annotation_id'
            for gff_file, annotation_id in queries:
                filtered_chunk = chunk[(chunk['gff_file'] == gff_file) & (chunk['annotation_id'] == annotation_id)]

                # Check if any rows pass the filter
                if not filtered_chunk.empty:
                    # Find the row with the longest sequence in the filtered chunk
                    max_length_index = filtered_chunk['dna_sequence'].str.len().idxmax()
                    max_length_sequence = filtered_chunk.loc[max_length_index, 'dna_sequence']

                    # Update the longest sequence and corresponding row if needed
                    if len(max_length_sequence) > max_sequence_length:
                        max_sequence_length = len(max_length_sequence)
                        max_sequence_row = filtered_chunk.loc[max_length_index]

                        # Check if the maximum length sequence is found
                        if max_sequence_length >= max_len_seq:
                            break

            # Check if the maximum length sequence is found
            if max_sequence_length >= max_len_seq:
                break

        # Get longest DNA sequence found
        if max_sequence_row is not None:
            longest_seq = max_sequence_row['dna_sequence']
            return longest_seq
        else:
            return None

    # Function to process each gene
    def process_gene(i, test_gene, igopan_notinfasta_len, igopan_all_panaroo_notinfasta, file_path):
        print(f"Processing accessory gene {i}: {test_gene}")

        # Filter to only get one of the genes of interest
        panaroo_seqs_interest_len = igopan_notinfasta_len[igopan_notinfasta_len.index == test_gene]
        panaroo_seqs_interest = igopan_all_panaroo_notinfasta[igopan_all_panaroo_notinfasta.index == test_gene]

        # Get length of longest sequence for that gene
        max_len_seq = panaroo_seqs_interest_len['Max group size nuc'][0]

        # Get a dataframe that has a list of gff file identifiers tied to the sequence id (sometimes more than one) found in that sample
        melted_df = panaroo_seqs_interest.melt(var_name='gff_file', value_name='seq_id')
        # Split rows with multiple values separated by ';' for the ones that have multiple sequences found in the same sample
        melted_df['seq_id'] = melted_df['seq_id'].str.split(';')
        melted_df = melted_df.explode('seq_id')
        melted_df = melted_df.dropna(subset=['seq_id'])
        melted_df = melted_df.reset_index(drop=True)
        queries = [(row.gff_file, row.seq_id) for row in melted_df.itertuples(index=False)]

        # Call function that gets the representative sequence
        rep_seq = find_rep_seq(file_path, queries, max_len_seq)

        return {'Gene': test_gene, 'rep_dna_sequence': rep_seq}

    if __name__ == '__main__':
        # Setup code and data initialization
        merge_presabs_accgenes = genefreq_meta_filt_acc.copy()

        # List that contains all the accessory COG IDs from earlier in this notebook
        shared_acccog_names = acccog_name_list

        # Obtain list of accessory genes in the fasta output file
        acc_genes_in_out = extract_sequence_identifiers(acc_fasta_path)

        # Obtain list of accessory genes that are NOT in the fasta output file
        acc_genes_notin_out = list(set(shared_acccog_names) - set(acc_genes_in_out))

        # Get gene_presence_absence panaroo output info of only the accessory genes not found in the panaroo fasta output file
        igopan_all_panaroo_notinfasta = igopan_all_panaroo.drop(['Non-unique Gene name', 'Annotation','No. isolates','No. sequences','Avg sequences per isolate','Genome Fragment','Order within Fragment','Accessory Fragment','Accessory Order with Fragment','QC','Min group size nuc','Avg group size nuc'], axis=1)
        igopan_all_panaroo_notinfasta.set_index(igopan_all_panaroo_notinfasta.columns[0], inplace=True)
        igopan_all_panaroo_notinfasta = igopan_all_panaroo_notinfasta[igopan_all_panaroo_notinfasta.index.isin(acc_genes_notin_out)]

        # Get dataframe with only max sequence length
        igopan_notinfasta_len = igopan_all_panaroo_notinfasta[['Max group size nuc']].copy()
        igopan_all_panaroo_notinfasta = igopan_all_panaroo_notinfasta.drop(['Max group size nuc'], axis=1)

        file_path = data_path + f'/gene_data.csv'

        # Parallelize process
        with ThreadPoolExecutor(max_workers=100) as executor:
            futures = [executor.submit(process_gene, i, acc_genes_notin_out[i], igopan_notinfasta_len, igopan_all_panaroo_notinfasta, file_path) for i in range(len(acc_genes_notin_out))]
            results = [future.result() for future in as_completed(futures)]

        # Convert results to DataFrame
        seqs_not_found = pd.DataFrame(results)

    ### Add the missing representative sequences to the exisisting acc gene fasta file

    # Function to append sequences to an existing FASTA file
    def append_to_fasta(df, fasta_path):
        # Create a list to hold SeqRecord objects
        records = []

        # Loop through the DataFrame rows
        for index, row in df.iterrows():
            # Create a SeqRecord object for each row
            seq_record = SeqRecord(
                Seq(row['rep_dna_sequence']),
                id=row['Gene'],
                description=""  # Empty description
            )
            records.append(seq_record)

        # Append the new records to the existing FASTA file
        with open(fasta_path, "a") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

    # Call the function to append the sequences to the existing FASTA file
    append_to_fasta(seqs_not_found, acc_fasta_path)

def get_pop_core_pres_abs_panaroo(data_path, out_path, core_lower):

    pres_abs_path = data_path+"/gene_presence_absence_roary.csv"
    sample_needed_path = data_path+'/needed_sample_names.txt'

    # Import panaroo output
    igopan_all_panaroo = pd.read_csv(pres_abs_path, low_memory=False)

    # Get list of panaroo output names in a list
    panpanaroo_ids_list =  list(igopan_all_panaroo["Gene"])

    # First round of filtering by fragment length >150bp OR general isolate frequency >10%

    tenp = (igopan_all_panaroo.shape[1]-14)/10 # Counting number of isolates (first 14 columns are metadata)

    panaroo_onefilt = igopan_all_panaroo[(igopan_all_panaroo['Avg group size nuc'] >= 150) | (igopan_all_panaroo['No. isolates'] >= tenp)]

    # Now make gene names the indeces
    panaroo_onefilt.set_index('Gene', inplace=True)

    # Drop all columns that are not an isolate
    panaroo_isol = panaroo_onefilt.iloc[:,13:]

    # Now replace NaN values for 0 and any other value for 1
    panaroo_isol[~panaroo_isol.isnull()] = 1
    panaroo_isol[panaroo_isol.isnull()] = 0

    # Switch rows to columns for ease
    panaroo_genefreq_matrix = panaroo_isol.transpose()
    panaroo_genefreq_matrix.index.name='Accession'

    # Get list of sample names in the subpopulation to analyse
    with open(sample_needed_path, 'r') as file:
        core_needed_list = file.read().splitlines()

    # Get only isolates in the list
    genefreq_mat_filt = panaroo_genefreq_matrix[panaroo_genefreq_matrix.index.isin(core_needed_list)]

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

    # Population all
    freq_cog = get_freq(genefreq_mat_filt)

    # Now filter to only keep cogs with frequency over a given value within the given sample subset
    core_cog = freq_cog.loc[(freq_cog['freq'] > core_lower)]

    # Get the number of core cogs per dataset
    corecog_num = core_cog.shape[0]

    # Get core cog names as lists
    corecog_name_list = list(core_cog["COG_name"])

    # Only get presence absence matrix for core genes
    genefreq_meta_filt_core = genefreq_mat_filt[genefreq_mat_filt.columns[genefreq_mat_filt.columns.isin(corecog_name_list)]]

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
            cog_id = record.description
            if cog_id in corecog_name_list:
                # Write the record to the output file
                record.id = cog_id
                record.description = ''
                SeqIO.write(record, outfile, "fasta")

    ### Find the rep sequences of the core genes not represented in the pan_genome_reference.fa Panaroo output

    # Function to write fasta files
    def fasta_out(rec):
        return(">"+str(rec.description)+"\n"+str(rec.seq))

    # Function to format the header for comparison
    def format(rec):
        return(rec.id)

    # Function to extract sequence identifiers from a fasta file
    def extract_sequence_identifiers(fasta_file_path):
        identifiers = []
        with open(fasta_file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    # Extract the sequence identifier (remove the ">")
                    identifier = line[1:]
                    identifiers.append(identifier)
        return identifiers

    # Function to find the representative sequence
    def find_rep_seq(file_path, queries, max_len_seq):
        chunk_iter = pd.read_csv(file_path, chunksize=1000)
        max_sequence_length = 0
        max_sequence_row = None

        # Loop through each chunk
        for chunk in chunk_iter:
            # Filter rows based on 'gff_file' and 'annotation_id'
            for gff_file, annotation_id in queries:
                filtered_chunk = chunk[(chunk['gff_file'] == gff_file) & (chunk['annotation_id'] == annotation_id)]

                # Check if any rows pass the filter
                if not filtered_chunk.empty:
                    # Find the row with the longest sequence in the filtered chunk
                    max_length_index = filtered_chunk['dna_sequence'].str.len().idxmax()
                    max_length_sequence = filtered_chunk.loc[max_length_index, 'dna_sequence']

                    # Update the longest sequence and corresponding row if needed
                    if len(max_length_sequence) > max_sequence_length:
                        max_sequence_length = len(max_length_sequence)
                        max_sequence_row = filtered_chunk.loc[max_length_index]

                        # Check if the maximum length sequence is found
                        if max_sequence_length >= max_len_seq:
                            break

            # Check if the maximum length sequence is found
            if max_sequence_length >= max_len_seq:
                break

        # Get longest DNA sequence found
        if max_sequence_row is not None:
            longest_seq = max_sequence_row['dna_sequence']
            return longest_seq
        else:
            return None

    # Function to process each gene
    def process_gene(i, test_gene, igopan_notinfasta_len, igopan_all_panaroo_notinfasta, file_path):
        print(f"Processing core gene {i}: {test_gene}")

        # Filter to only get one of the genes of interest
        panaroo_seqs_interest_len = igopan_notinfasta_len[igopan_notinfasta_len.index == test_gene]
        panaroo_seqs_interest = igopan_all_panaroo_notinfasta[igopan_all_panaroo_notinfasta.index == test_gene]

        # Get length of longest sequence for that gene
        max_len_seq = panaroo_seqs_interest_len['Max group size nuc'][0]

        # Get a dataframe that has a list of gff file identifiers tied to the sequence id (sometimes more than one) found in that sample
        melted_df = panaroo_seqs_interest.melt(var_name='gff_file', value_name='seq_id')
        # Split rows with multiple values separated by ';' for the ones that have multiple sequences found in the same sample
        melted_df['seq_id'] = melted_df['seq_id'].str.split(';')
        melted_df = melted_df.explode('seq_id')
        melted_df = melted_df.dropna(subset=['seq_id'])
        melted_df = melted_df.reset_index(drop=True)
        queries = [(row.gff_file, row.seq_id) for row in melted_df.itertuples(index=False)]

        # Call function that gets the representative sequence
        rep_seq = find_rep_seq(file_path, queries, max_len_seq)

        return {'Gene': test_gene, 'rep_dna_sequence': rep_seq}

    if __name__ == '__main__':
        # Setup code and data initialization
        merge_presabs_coregenes = genefreq_meta_filt_core.copy()

        # List that contains all the core COG IDs from earlier in this notebook
        shared_corecog_names = corecog_name_list

        # Obtain list of core genes in the fasta output file
        core_genes_in_out = extract_sequence_identifiers(core_fasta_path)

        # Obtain list of core genes that are NOT in the fasta output file
        core_genes_notin_out = list(set(shared_corecog_names) - set(core_genes_in_out))

        # Get gene_presence_absence panaroo output info of only the core genes not found in the panaroo fasta output file
        igopan_all_panaroo_notinfasta = igopan_all_panaroo.drop(['Non-unique Gene name', 'Annotation','No. isolates','No. sequences','Avg sequences per isolate','Genome Fragment','Order within Fragment','Accessory Fragment','Accessory Order with Fragment','QC','Min group size nuc','Avg group size nuc'], axis=1)
        igopan_all_panaroo_notinfasta.set_index(igopan_all_panaroo_notinfasta.columns[0], inplace=True)
        igopan_all_panaroo_notinfasta = igopan_all_panaroo_notinfasta[igopan_all_panaroo_notinfasta.index.isin(core_genes_notin_out)]

        # Get dataframe with only max sequence length
        igopan_notinfasta_len = igopan_all_panaroo_notinfasta[['Max group size nuc']].copy()
        igopan_all_panaroo_notinfasta = igopan_all_panaroo_notinfasta.drop(['Max group size nuc'], axis=1)

        file_path = data_path + f'/gene_data.csv'

        # Use ThreadPoolExecutor to parallelize the process
        with ThreadPoolExecutor(max_workers=100) as executor:
            futures = [executor.submit(process_gene, i, core_genes_notin_out[i], igopan_notinfasta_len, igopan_all_panaroo_notinfasta, file_path) for i in range(len(core_genes_notin_out))]
            results = [future.result() for future in as_completed(futures)]

        # Convert results to DataFrame
        seqs_not_found = pd.DataFrame(results)

    ### Add the missing representative sequences to the exisisting core gene fasta file

    # Function to append sequences to an existing FASTA file
    def append_to_fasta(df, fasta_path):
        # Create a list to hold SeqRecord objects
        records = []

        # Loop through the DataFrame rows
        for index, row in df.iterrows():
            # Create a SeqRecord object for each row
            seq_record = SeqRecord(
                Seq(row['rep_dna_sequence']),
                id=row['Gene'],
                description=""  # Empty description
            )
            records.append(seq_record)

        # Append the new records to the existing FASTA file
        with open(fasta_path, "a") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

    # Call the function to append the sequences to the existing FASTA file
    append_to_fasta(seqs_not_found, core_fasta_path)

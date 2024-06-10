#!/usr/bin/env python
# coding: utf-8

# EggNOG database wrapper

# Since EggNOG-mapper is down (rip), and I do not want for my pipeline to depend on a website that might be down, I am trying to build my own wrapper of the EggNOG database.
#
# What I ended up doing was to use the eggnog 5.0 database website (http://eggnog5.embl.de/#/app/home) which provides a protein sequence query option to functionally annotate my sequences, since they will always be in the order of 2,000-3,000 COGs, which should not take more than an hour to run. EggNOG-mapper uses a slightly different method that relies on fast alignment to seed orthologs using diamond, but I tried to get that running and many issues came up. Also it has the drawback of having to download the databases which take a lot of space. It does seem like the functional group distributions are comparable with only slight variations.
#
# I do think these annotations might be more accurate, since I am using the most abundant functional group annotated for the top 5 ortholog groups hits with the database. In the future I could refine it to prioritize hits with Streptococcus, other paramaters, etc. For now, this seemed like a good basic annotation.

# Hold this for now in case it doesn't work: "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:125.0) Gecko/20100101 Firefox/125.0"

# ### Import necessary packages

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import pandas as pd
import re
import itertools
import time
import os

# ### EggNOG annotation through the website database querying

def get_functional_groups(out_path):

    #fasta_path = out_path+'/accessory_rep_seqs.fasta'
    fasta_path = out_path+'/accessory_rep_seqs.fasta'

    # Function to submit sequences to the eggnog database website

    def submit_sequence(sequence, max_retries=3, delay=1):
        url = "http://eggnogapi5.embl.de/fast_webscan"
        headers = {
            "User-Agent": "Mozilla/5.0",
            "Accept": "application/json, text/plain, */*",
            "Accept-Language": "en-US,en;q=0.5",
            "Accept-Encoding": "gzip, deflate",
            "Content-Type": "application/json;charset=utf-8",
            "Origin": "http://eggnog5.embl.de",
            "Connection": "keep-alive",
            "Referer": "http://eggnog5.embl.de/",
            "DNT": "1",
            "Sec-GPC": "1"
        }
        data = {
            "fasta": sequence,
            "nqueries": 1,
            "db": "euk"
        }

        retries = 0
        while retries < max_retries:
            response = requests.post(url, headers=headers, json=data)
            if response.status_code == 200:
                return response.json()
            else:
                print(f"Failed to submit sequence. Retrying ({retries+1}/{max_retries})...")
                retries += 1
                time.sleep(delay)

        print("Max retries reached. Failed to submit sequence.")
        return None

    # Function to translate DNA sequences to protein sequences

    def translate_dna_to_protein(dna_sequence):
        dna_seq = Seq(dna_sequence)
        protein_seq = dna_seq.translate()
        # Remove the asterisk (*) at the end, which symbolizes a stop codon but confuses the database
        if protein_seq.endswith("*"):
            protein_seq = protein_seq[:-1]
        return str(protein_seq)

    # Extract most common functional groups from the top 5 hits in the eggnog database website query

    def extract_functional_group(sequence_result):
        seq_matches = sequence_result.get('seq_matches', [])
        funcat_counts = {}
        max_funcat = None
        max_count = 0

        for match in seq_matches:
            hits = match.get('hit', {}).get('matches', [])
            # Consider only the top 5 hits
            for hit in hits[:5]:
                funcat = hit.get('funcat')
                funcat_counts[funcat] = funcat_counts.get(funcat, 0) + 1
                if funcat_counts[funcat] > max_count:
                    max_count = funcat_counts[funcat]
                    max_funcat = funcat

        return max_funcat

    # Make dataframe for results
    eggnog_functional = pd.DataFrame(columns=['cog','prot_sequence','functional_group'])
    counter = -1

    sequences = SeqIO.parse(fasta_path, "fasta")

    start_time = time.time()  # Start timing

    # Iterate over sequences, translate, and extract functional group
    for sequence in sequences:
        counter = counter + 1

        translated_sequence = translate_dna_to_protein(str(sequence.seq))
        sequence_result = submit_sequence(translated_sequence)
        functional_group = extract_functional_group(sequence_result)
        if functional_group:
            funct_group = list(functional_group)[0]
        else:
            funct_group = 'NH'

        eggnog_functional.at[counter, 'cog'] = sequence.id
        eggnog_functional.at[counter, 'functional_group'] = funct_group
        eggnog_functional.at[counter, 'prot_sequence'] = translated_sequence

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Annotations completed in "+str(elapsed_time)+" seconds")

    # Create linkage folder to put the results there
    # Specify the directory path you want to create
    eggnog_path = out_path+"/eggnog"

    # Check if the directory already exists
    if not os.path.exists(eggnog_path):
        # Create the directory
        os.makedirs(eggnog_path)

    eggnog_functional.to_csv(eggnog_path+'/acc_cog_eggnog_annotations.csv', index=False)

    # Now to summarize the COGs that are part of each eggnog functional group

    # Only get COG names with exact match

    # A
    eggnog_acc_A = eggnog_functional[eggnog_functional["functional_group"] == 'A']
    eggnog_hit_cogs_A =  list(eggnog_acc_A["cog"])

    # C
    eggnog_acc_C = eggnog_functional[eggnog_functional["functional_group"] == 'C']
    eggnog_hit_cogs_C =  list(eggnog_acc_C["cog"])

    # D
    eggnog_acc_D = eggnog_functional[eggnog_functional["functional_group"] == 'D']
    eggnog_hit_cogs_D =  list(eggnog_acc_D["cog"])

    # E
    eggnog_acc_E = eggnog_functional[eggnog_functional["functional_group"] == 'E']
    eggnog_hit_cogs_E =  list(eggnog_acc_E["cog"])

    # F
    eggnog_acc_F = eggnog_functional[eggnog_functional["functional_group"] == 'F']
    eggnog_hit_cogs_F =  list(eggnog_acc_F["cog"])

    # G
    eggnog_acc_G = eggnog_functional[eggnog_functional["functional_group"] == 'G']
    eggnog_hit_cogs_G =  list(eggnog_acc_G["cog"])

    # H
    eggnog_acc_H = eggnog_functional[eggnog_functional["functional_group"] == 'H']
    eggnog_hit_cogs_H =  list(eggnog_acc_H["cog"])

    # I
    eggnog_acc_I = eggnog_functional[eggnog_functional["functional_group"] == 'I']
    eggnog_hit_cogs_I =  list(eggnog_acc_I["cog"])

    # J
    eggnog_acc_J = eggnog_functional[eggnog_functional["functional_group"] == 'J']
    eggnog_hit_cogs_J =  list(eggnog_acc_J["cog"])

    # K
    eggnog_acc_K = eggnog_functional[eggnog_functional["functional_group"] == 'K']
    eggnog_hit_cogs_K =  list(eggnog_acc_K["cog"])

    # L
    eggnog_acc_L = eggnog_functional[eggnog_functional["functional_group"] == 'L']
    eggnog_hit_cogs_L =  list(eggnog_acc_L["cog"])

    # M
    eggnog_acc_M = eggnog_functional[eggnog_functional["functional_group"] == 'M']
    eggnog_hit_cogs_M =  list(eggnog_acc_M["cog"])

    # N
    eggnog_acc_N = eggnog_functional[eggnog_functional["functional_group"] == 'N']
    eggnog_hit_cogs_N =  list(eggnog_acc_N["cog"])

    # O
    eggnog_acc_O = eggnog_functional[eggnog_functional["functional_group"] == 'O']
    eggnog_hit_cogs_O =  list(eggnog_acc_O["cog"])

    # P
    eggnog_acc_P = eggnog_functional[eggnog_functional["functional_group"] == 'P']
    eggnog_hit_cogs_P =  list(eggnog_acc_P["cog"])

    # Q
    eggnog_acc_Q = eggnog_functional[eggnog_functional["functional_group"] == 'Q']
    eggnog_hit_cogs_Q =  list(eggnog_acc_Q["cog"])

    # T
    eggnog_acc_T = eggnog_functional[eggnog_functional["functional_group"] == 'T']
    eggnog_hit_cogs_T =  list(eggnog_acc_T["cog"])

    # U
    eggnog_acc_U = eggnog_functional[eggnog_functional["functional_group"] == 'U']
    eggnog_hit_cogs_U =  list(eggnog_acc_U["cog"])

    # V
    eggnog_acc_V = eggnog_functional[eggnog_functional["functional_group"] == 'V']
    eggnog_hit_cogs_V =  list(eggnog_acc_V["cog"])

    # S
    eggnog_acc_S = eggnog_functional[eggnog_functional["functional_group"] == 'S']
    eggnog_hit_cogs_S =  list(eggnog_acc_S["cog"])

    # NH
    eggnog_acc_NH = eggnog_functional[eggnog_functional["functional_group"] == 'NH']
    eggnog_hit_cogs_NH =  list(eggnog_acc_NH["cog"])

    # Mixed
    eggnog_hit_cogs_mixed = []
    for i in range(len(eggnog_functional)):
        if eggnog_functional.iloc[i]["functional_group"] != 'NH':
            if eggnog_functional.functional_group.str.len()[i] > 1:
                cog = eggnog_functional.iloc[i]["cog"]
                eggnog_hit_cogs_mixed.append(cog)


    # Create a dataframe per cog functional group
    # All
    acc_gene_names = list(eggnog_functional["cog"])
    all_cogs = pd.DataFrame({'all_acc_COGs': acc_gene_names})

    # Eggnog
    eggnog_A = pd.DataFrame({'eggnog_A': sorted(eggnog_hit_cogs_A)})
    eggnog_C = pd.DataFrame({'eggnog_C': sorted(eggnog_hit_cogs_C)})
    eggnog_D = pd.DataFrame({'eggnog_D': sorted(eggnog_hit_cogs_D)})
    eggnog_E = pd.DataFrame({'eggnog_E': sorted(eggnog_hit_cogs_E)})
    eggnog_F = pd.DataFrame({'eggnog_F': sorted(eggnog_hit_cogs_F)})
    eggnog_G = pd.DataFrame({'eggnog_G': sorted(eggnog_hit_cogs_G)})
    eggnog_H = pd.DataFrame({'eggnog_H': sorted(eggnog_hit_cogs_H)})
    eggnog_I = pd.DataFrame({'eggnog_I': sorted(eggnog_hit_cogs_I)})
    eggnog_J = pd.DataFrame({'eggnog_J': sorted(eggnog_hit_cogs_J)})
    eggnog_K = pd.DataFrame({'eggnog_K': sorted(eggnog_hit_cogs_K)})
    eggnog_L = pd.DataFrame({'eggnog_L': sorted(eggnog_hit_cogs_L)})
    eggnog_M = pd.DataFrame({'eggnog_M': sorted(eggnog_hit_cogs_M)})
    eggnog_N = pd.DataFrame({'eggnog_N': sorted(eggnog_hit_cogs_N)})
    eggnog_O = pd.DataFrame({'eggnog_O': sorted(eggnog_hit_cogs_O)})
    eggnog_P = pd.DataFrame({'eggnog_P': sorted(eggnog_hit_cogs_P)})
    eggnog_Q = pd.DataFrame({'eggnog_Q': sorted(eggnog_hit_cogs_Q)})
    eggnog_T = pd.DataFrame({'eggnog_T': sorted(eggnog_hit_cogs_T)})
    eggnog_U = pd.DataFrame({'eggnog_U': sorted(eggnog_hit_cogs_U)})
    eggnog_V = pd.DataFrame({'eggnog_V': sorted(eggnog_hit_cogs_V)})
    eggnog_S = pd.DataFrame({'eggnog_S': sorted(eggnog_hit_cogs_S)})
    eggnog_NH = pd.DataFrame({'eggnog_NH': sorted(eggnog_hit_cogs_NH)})
    eggnog_mixed = pd.DataFrame({'eggnog_mixed': sorted(eggnog_hit_cogs_mixed)})

    # Concatenate dataframes
    functionalCOG_namelist = pd.concat([all_cogs,eggnog_A,eggnog_C,eggnog_D,eggnog_E,eggnog_F,eggnog_G,eggnog_H,eggnog_I,eggnog_J,eggnog_K,eggnog_L,eggnog_M,eggnog_N,eggnog_O,eggnog_P,eggnog_Q,eggnog_T,eggnog_U,eggnog_V,eggnog_S,eggnog_NH,eggnog_mixed], axis=1)

    functionalCOG_namelist.to_csv(eggnog_path+'/eggnog_group_cog_names.csv', index=False)

    # Export fasta with protein sequences for the accessory genes
    protein_fasta_path = out_path+'/accessory_rep_protein_seqs.fasta'
    with open(protein_fasta_path, "w") as fasta_file:
        for index, row in eggnog_functional.iterrows():
            cog = row['cog']
            sequence = row['prot_sequence']
            fasta_file.write(f">{cog}\n{sequence}\n")

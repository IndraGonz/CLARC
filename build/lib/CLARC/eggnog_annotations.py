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

import os
import time
import random
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq

# Starting from CLARC v1.1.3
# Revision list:
#- Parallelized submission process to the EggNOG API, choosing the optimal amount of parallel jobs from the API response time
#- Made sure to put latent time between the jobs as to not exceed the request limits (this caps the efficiency but it's good practice)
#- Made the commands to summarize the COG functional category results more efficient

# ### EggNOG annotation through the website database querying

def get_functional_groups(out_path, max_cores):

    """Main function to process sequences and store EggNOG annotations."""

    session = requests.Session()
    max_cores = max(1, max_cores // 2)  # Ensure at least 1 core is used

    def submit_sequence(sequence, max_retries=3):
        """Submits a protein sequence to EggNOG API with retry logic."""
        url = "http://eggnogapi5.embl.de/fast_webscan"
        headers = {
            "User-Agent": "Mozilla/5.0",
            "Accept": "application/json, text/plain, */*",
            "Content-Type": "application/json;charset=utf-8",
        }
        data = {"fasta": sequence, "nqueries": 1, "db": "bact"}

        for attempt in range(max_retries):
            try:
                response = session.post(url, headers=headers, json=data)
                if response.status_code == 200:
                    return response.json()
                elif response.status_code == 429:  # Too Many Requests
                    wait_time = 2 ** attempt  # Exponential backoff
                    print(f"Rate limit hit! Waiting {wait_time} sec...")
                    time.sleep(wait_time)
                elif response.status_code in [500, 503]:
                    wait_time = random.uniform(1, 3)
                    print(f"Server error {response.status_code}. Retrying in {wait_time:.1f} sec...")
                    time.sleep(wait_time)
                else:
                    print(f"Failed request (status {response.status_code}). Retrying...")

            except requests.exceptions.RequestException as e:
                print(f"Request error: {e}")

            time.sleep(random.uniform(0.05, 0.1))

        print("Max retries reached; failed to submit sequence.")
        return None

    def translate_dna_to_protein(dna_sequence):
        """Translates DNA to a protein sequence, removing stop codon (*)"""
        return str(Seq(dna_sequence).translate()).rstrip("*")

    def extract_functional_group(sequence_result):
        """Extracts the most common functional group from the EggNOG API response."""
        if not sequence_result:
            return 'NH'

        funcat_counts = {}
        for match in sequence_result.get('seq_matches', []):
            for hit in match.get('hit', {}).get('matches', [])[:5]:
                funcat = hit.get('funcat')
                if funcat:
                    funcat_counts[funcat] = funcat_counts.get(funcat, 0) + 1

        return max(funcat_counts, key=funcat_counts.get, default='NH')

    def process_sequence(seq):
        """Processes a single sequence: translates DNA, submits to API, extracts functional group."""
        translated_sequence = translate_dna_to_protein(str(seq.seq))
        sequence_result = submit_sequence(translated_sequence)
        functional_group = extract_functional_group(sequence_result)

        return {
            'cog': seq.id,
            'prot_sequence': translated_sequence,
            'functional_group': functional_group
        }

    def measure_api_speed():
        """Tests API response time to determine optimal concurrency."""
        test_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        translated_seq = translate_dna_to_protein(test_seq)

        start_time = time.time()
        response = submit_sequence(translated_seq)
        elapsed_time = time.time() - start_time

        return elapsed_time if response else None

    def auto_tune_n_jobs():
        """Dynamically adjusts parallel job count based on API speed."""
        times = [measure_api_speed() for _ in range(3)]
        avg_time = sum(times) / len(times) if all(times) else None

        if avg_time is None:
            return 2

        if avg_time < 1:
            return min(20, max_cores)
        elif avg_time < 3:
            return min(12, max_cores)
        elif avg_time < 5:
            return min(6, max_cores)
        else:
            return 2

    ### Start Processing ###
    fasta_path = os.path.join(out_path, 'accessory_rep_seqs.fasta')

    if not os.path.exists(fasta_path):
        print(f"Error: Input file '{fasta_path}' not found.")
        return

    sequences = list(SeqIO.parse(fasta_path, "fasta"))

    if not sequences:
        print("Error: No sequences found in the input file.")
        return

    start_time = time.time()

    optimal_n_jobs = auto_tune_n_jobs()

    results = []
    with ThreadPoolExecutor(max_workers=optimal_n_jobs) as executor:
        future_to_seq = {executor.submit(process_sequence, seq): seq.id for seq in sequences}
        for future in as_completed(future_to_seq):
            results.append(future.result())

    eggnog_functional = pd.DataFrame(results)

    eggnog_path = os.path.join(out_path, "eggnog")
    os.makedirs(eggnog_path, exist_ok=True)
    eggnog_functional.to_csv(os.path.join(eggnog_path, 'acc_cog_eggnog_annotations.csv'), index=False)

    unique_groups = set("".join(eggnog_functional["functional_group"].dropna().unique()))
    unique_groups.add("NH")
    unique_groups.add("eggnog_mixed")

    grouped_cogs = {f"eggnog_{g}": [] for g in sorted(unique_groups)}
    grouped_cogs["all_acc_COGs"] = []
    grouped_cogs["eggnog_mixed"] = []

    for _, row in eggnog_functional.iterrows():
        cog = row["cog"]
        func_group = row["functional_group"]

        grouped_cogs["all_acc_COGs"].append(cog)

        if len(func_group) > 1 and func_group not in ["NH", "NA"]:
            grouped_cogs["eggnog_mixed"].append(cog)
        else:
            grouped_cogs.setdefault(f"eggnog_{func_group}", []).append(cog)

    functionalCOG_namelist = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in grouped_cogs.items()]))
    functionalCOG_namelist.to_csv(os.path.join(eggnog_path, 'eggnog_group_cog_names.csv'), index=False)

    protein_fasta_path = os.path.join(out_path, 'accessory_rep_protein_seqs.fasta')
    with open(protein_fasta_path, "w") as fasta_file:
        for row in results:
            fasta_file.write(f">{row['cog']}\n{row['prot_sequence']}\n")

    print(f"EggNOG annotations completed in {time.time() - start_time:.2f} seconds")

#! /bin/bash

# Assign the first argument to out_path
out_path=${1}

# Define the directory path for BLAST output
blast_path="${out_path}/acc_blastn"

# Check if the directory exists, create if not
if [ ! -d "$blast_path" ]; then
    mkdir -p "$blast_path"
fi

# Copy file to be the database
cp "${out_path}/accessory_rep_seqs.fasta" "${blast_path}/accessory_rep_seqs_db.fasta"

# Set file/folder paths
acc_cogs_db="${blast_path}/accessory_rep_seqs_db.fasta"
acc_cogs="${out_path}/accessory_rep_seqs.fasta"

# Make databases
makeblastdb -in "${acc_cogs_db}" -dbtype nucl

# Run BLASTN search
blastn -task blastn -query "${acc_cogs}" -db "${acc_cogs_db}" -evalue 1e-20 -num_threads 4 -out "${blast_path}/blastn_acccogs_allvall.tsv" -outfmt 6

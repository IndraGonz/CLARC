#!/usr/bin/env python

# Author: Indra González Ojeda
# CLARC: Pipeline that cleans the COG definitions of a specific population within a multi-population pangenome analysis

import argparse
import sys
import os
import shutil
from .filtering_acc_core import get_pop_acc_pres_abs, get_pop_core_pres_abs
from .filtering_acc_core_panaroo import get_pop_acc_pres_abs_panaroo, get_pop_core_pres_abs_panaroo
from .get_linkage_matrix import get_linkage_matrices
from .eggnog_annotations import get_functional_groups
from .clarc_condense import clarc_cleaning
from .clarc_merge import merge_clarc_results
import subprocess
from .version import __version__  # Import the version

# Run stuff

def main():
    parser = argparse.ArgumentParser(description="CLARC pipeline to clean COG definitions for individual populations within a general pangenome analysis")
    parser.add_argument("--input_dir", help="Directory that contains input data")
    parser.add_argument("--output_dir", default="output", help="Output directory (results saved here)")
    parser.add_argument("--acc_upper", default=0.95, type=float, help="Upper bound for accessory gene filtering")
    parser.add_argument("--acc_lower", default=0.05, type=float, help="Lower bound for accessory gene filtering")
    parser.add_argument("--core_lower", default=0.95, type=float, help="Lower bound for core gene filtering")
    parser.add_argument("--ci","--clarc-identity", default=95, type=float, help="BLASTN identity threshold CLARC uses as constraint to identity same gene clusters. Number from 0-100, default is 95%")
    parser.add_argument("--filter-only", action='store_true', help="If given, only run the filtering steps")
    parser.add_argument("--panaroo", action='store_true', help="If given, the input data will be from Panaroo and it will filter accordingly. Remember to provide the 'gene_data.csv' input")
    parser.add_argument("--options", action='store_true', help="Show available options")
    parser.add_argument("--version", action='store_true', help="Show the version of the CLARC tool")
    parser.add_argument("--dif", "--delete-intermediate-files", action='store_true', help="Delete intermediate files")
    parser.add_argument('--merge', nargs='+', help='Paths to CLARC results to merge (folders where clarc was run that include the original data and clarc output subfolders)')

    args = parser.parse_args()

    # Error if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Display options if --options is called
    if args.options:
        parser.print_help()
        sys.exit(0)

    # Display version if --version is called
    if args.version:
        print(f"CLARC version {__version__}")
        sys.exit(0)

    ## Set input parameters and PATHs ####
    input_dir = args.input_dir
    output_dir = args.output_dir
    acc_upper = args.acc_upper
    acc_lower = args.acc_lower
    core_lower = args.core_lower
    filter_only = args.filter_only
    panaroo = args.panaroo
    merge = args.merge
    ci = args.ci

    if merge:

        print("Merging results from different CLARC runs to obtain final cluster list...")

        if args.panaroo:
            panaroo_true = 1
        else:
            panaroo_true = 0

        merge_clarc_results(merge, out_path=args.output_dir, panaroo_true=panaroo_true)

    else:

        print("Filtering data to identify accessory and core genes...")

        if args.panaroo:
            # Filter Panaroo data to create presence absence matrices for core and accessory genes
            panaroo_true = 1
            get_pop_acc_pres_abs_panaroo(input_dir, output_dir, acc_upper, acc_lower)
            get_pop_core_pres_abs_panaroo(input_dir, output_dir, core_lower)

        else:
            panaroo_true = 0
            # Filter Roary data to create presence absence matrices for core and accessory genes
            get_pop_acc_pres_abs(input_dir, output_dir, acc_upper, acc_lower)
            get_pop_core_pres_abs(input_dir, output_dir, core_lower)

        print("Data filtered, presence absence matrices created for subpopulation of samples.")

        # If --filter-only is given, stop here
        if filter_only:
            print("CLARC finished running on 'filter only' mode")
            return

        print("Calculating the linkage matrices...")

        ## Get linkage matrices for CLARC analysis
        get_linkage_matrices(output_dir)

        print("Linkage matrices generated for the subpopulation accessory genes.")

        ## Perform blastn all v all comparison for accessory gene rep sequences
        subprocess.run(['bash', 'acccog_blastn.sh', output_dir])

        print("All vs. all nucleotide BLAST performed for the subpopulation accessory genes.")

        print("Performing EggNOG functional annotations...")

        ## Perform eggnog functional annotation
        get_functional_groups(output_dir)

        print("Eggnog functional annotation of accessory COGs is complete.")

        ## Perform CLARC cleaning of COGs
        clarc_cleaning(input_dir, output_dir, panaroo_true, acc_upper, acc_lower, core_lower, ci)

        print("CLARC finished re-defining COGs. Have fun with the results!")

        if args.dif:
            # List of paths to delete intermediate files
            paths_to_delete = [
                output_dir + '/linkage',
                output_dir + '/eggnog',
                output_dir + '/acc_blastn',
                output_dir + '/accessory_rep_protein_seqs.fasta',
                output_dir + '/clarc_results/accessory_cluster_cogs.txt',
                output_dir + '/clarc_results/accessory_cluster_summary.csv',
                output_dir + '/clarc_results/core_cluster_cogs.txt',
                output_dir + '/clarc_results/core_cluster_summary.csv'
            ]

            # Iterate through each path and delete accordingly
            for path in paths_to_delete:
                try:
                    if os.path.isfile(path):
                        os.remove(path)
                    elif os.path.isdir(path):
                        shutil.rmtree(path)
                except Exception as e:
                    print(f"Error deleting {path}: {e}")


if __name__ == "__main__":
    main()

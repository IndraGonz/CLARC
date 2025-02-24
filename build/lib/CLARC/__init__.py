#!/usr/bin/env python

# Author: Indra González Ojeda
# CLARC: Pipeline that cleans the COG definitions of a specific population within a multi-population pangenome analysis

import argparse
import sys
import os
import shutil
from .filtering_acc_core import get_pop_acc_pres_abs, get_pop_core_pres_abs
from .filtering_acc_core_panaroo import get_pop_acc_pres_abs_panaroo, get_pop_core_pres_abs_panaroo
from .filtering_acc_core_ppanggo import get_pop_acc_pres_abs_ppanggo, get_pop_core_pres_abs_ppanggo
from .get_linkage_matrix import get_linkage_matrices
from .eggnog_annotations import get_functional_groups
from .clarc_condense import clarc_cleaning
from .clarc_merge import merge_clarc_results
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from .version import __version__  # Import the version

# Run stuff

def main():
    parser = argparse.ArgumentParser(description="CLARC pipeline to clean COG definitions for individual populations within a general pangenome analysis")
    parser.add_argument("--input_dir", help="Directory that contains input data")
    parser.add_argument("--output_dir", default="clarc_output", help="Output directory (results saved here)")
    parser.add_argument("--acc_upper", default=0.95, type=float, help="Upper bound for accessory gene filtering")
    parser.add_argument("--acc_lower", default=0.05, type=float, help="Lower bound for accessory gene filtering")
    parser.add_argument("--core_lower", default=0.95, type=float, help="Lower bound for core gene filtering")
    parser.add_argument("--linkage_cut", default=0, type=float, help="Linkage threshold for co-occurrence constraint. Default is 0 which means COG pairs that fully exclude each other")
    parser.add_argument("--connection_cut", default=1, type=float, help="Threshold for connectivity of clusters to condense. Default is 1, only fully connected clusters are condensed")
    parser.add_argument("--ci","--clarc-identity", default=95, type=float, help="BLASTN identity threshold CLARC uses as constraint to identity same gene clusters. Number from 0-100, default is 95 percent")
    parser.add_argument("--filter-only", action='store_true', help="If given, only run the filtering steps")
    parser.add_argument("--panaroo", action='store_true', help="If given, the input data will be from Panaroo and it will filter accordingly. Remember to provide the 'gene_data.csv' input")
    parser.add_argument("--ppanggo", action='store_true', help="If given, the input data will be from PPanGGolin and it will filter accordingly. Remember to provide the appropiate inputs")
    parser.add_argument("--options", action='store_true', help="Show available options")
    parser.add_argument("--version", action='store_true', help="Show the version of the CLARC tool")
    parser.add_argument("--max_cores", type=int, default=os.cpu_count(), help="Maximum number of CPU cores to use for parallel tasks. Default is all available cores")
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
    ppanggo = args.ppanggo
    merge = args.merge
    ci = args.ci
    linkage_cut = args.linkage_cut
    connection_cut = args.connection_cut
    max_cores = args.max_cores

    # Ensure the output directory is handled correctly
    if os.path.exists(output_dir):
        print(f"Warning: Output directory '{output_dir}' already exists. Files may be overwritten")
    else:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")

    ## Determine number of CPU cores to use
    max_cores = min(args.max_cores, os.cpu_count())  # Ensure we don't exceed available CPUs
    print(f"Using {max_cores} cores for parallel processing")

    ## Set paths to intermediate files to check if the step has been previously run
    filter_path = os.path.join(output_dir, "population_accessory_presence_absence.csv")
    linkage_path = os.path.join(output_dir, "linkage", "acc_linkage_co-occur.csv")
    eggnog_path = os.path.join(output_dir, "eggnog", "acc_cog_eggnog_annotations.csv")
    blast_path = os.path.join(output_dir, "acc_blastn", "blastn_acccogs_allvall.tsv")


    if args.panaroo:
        panaroo_true = 1
    else:
        panaroo_true = 0

    if args.ppanggo:
        ppanggo_true = 1
    else:
        ppanggo_true = 0

    if merge:

        print("Merging results from different CLARC runs to obtain final cluster list...")

        merge_clarc_results(merge, out_path=args.output_dir, panaroo_true=panaroo_true)

    else:

        def filter():
            if not os.path.exists(filter_path):

                print("Filtering data to identify accessory and core genes...")

                if panaroo_true == 1:

                    get_pop_acc_pres_abs_panaroo(input_dir, output_dir, acc_upper, acc_lower)
                    get_pop_core_pres_abs_panaroo(input_dir, output_dir, core_lower)

                elif ppanggo_true == 1:

                    get_pop_acc_pres_abs_ppanggo(input_dir, output_dir, acc_upper, acc_lower)
                    get_pop_core_pres_abs_ppanggo(input_dir, output_dir, core_lower)

                else:

                    # Filter Roary data to create presence absence matrices for core and accessory genes
                    get_pop_acc_pres_abs(input_dir, output_dir, acc_upper, acc_lower)
                    get_pop_core_pres_abs(input_dir, output_dir, core_lower)

                    print("Data filtered, presence absence matrices created for subpopulation of samples")
            else:
                print("Skipping filtering step, output already exists")

        # If --filter-only is given, stop here
        if filter_only:
            print("CLARC finished running on 'filter only' mode")
            return

        def linkage():
            if not os.path.exists(linkage_path):

                print("Calculating the linkage matrices...")

                ## Get linkage matrices for CLARC analysis
                get_linkage_matrices(output_dir, max_cores)

                print("Linkage matrices generated for the subpopulation accessory genes")

            else:
                print("Skipping linkage matrix calculation step, output already exists")

        def blastn():
            if not os.path.exists(blast_path):

                ## Perform blastn all v all comparison for accessory gene rep sequences
                subprocess.run(['bash', 'acccog_blastn.sh', output_dir])

                print("All vs. all nucleotide BLAST performed for the subpopulation accessory genes")

            else:
                print("Skipping all vs all blast step, output already exists")

        def eggnog():
            if not os.path.exists(eggnog_path):

                print("Performing EggNOG functional annotations...")

                ## Perform eggnog functional annotation
                get_functional_groups(output_dir, max_cores)

                print("Eggnog functional annotation of accessory COGs is complete")


            else:
                print("Skipping EggNOG functional annotation step, output already exists")

        ## Call functions to perform the CLARC preparation steps (and check if these steps have been previously run)
        filter()

        # Run linkage, blastn, and eggnog in parallel
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {
                executor.submit(linkage): "Linkage",
                executor.submit(blastn): "BLASTN",
                executor.submit(eggnog): "EggNOG"
            }

            for future in as_completed(futures):
                task_name = futures[future]
                try:
                    future.result()
                    print(f"{task_name} step completed successfully.")
                except Exception as e:
                    print(f"Error in {task_name} step: {e}")
                    sys.exit(1)

        ## Perform CLARC cleaning of COGs
        clarc_cleaning(input_dir, output_dir, panaroo_true, ppanggo_true, acc_upper, acc_lower, core_lower, ci, linkage_cut, connection_cut)

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

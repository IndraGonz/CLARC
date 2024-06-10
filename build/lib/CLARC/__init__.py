#!/usr/bin/env python

# Author: Indra González Ojeda
# CLARC: Pipeline that cleans gene definitions within a population in a multi-population pangenome analysis

import argparse
from .filtering_acc_core import get_pop_acc_pres_abs, get_pop_core_pres_abs
from .get_linkage_matrix import get_linkage_matrices
from .eggnog_annotations import get_functional_groups
from .clarc_condense import clarc_cleaning
import pandas as pd
import subprocess

# Run stuff

def main():
    parser = argparse.ArgumentParser(description="CLARC pipeline to clean COG definitions for individual populations within a general pangenome analysis")
    parser.add_argument("--input_dir", help="Directory that contains input data")
    parser.add_argument("--output_dir", default="output", help="Output directory (results saved here)")
    parser.add_argument("--acc_upper", default=0.95, type=float, help="Upper bound for accessory gene filtering")
    parser.add_argument("--acc_lower", default=0.05, type=float, help="Lower bound for accessory gene filtering")
    parser.add_argument("--core_lower", default=0.95, type=float, help="Lower bound for core gene filtering")
    args = parser.parse_args()

    ## Set input parameters and PATHs ####
    input_dir = args.input_dir
    output_dir = args.output_dir
    acc_upper = args.acc_upper
    acc_lower = args.acc_lower
    core_lower = args.core_lower

    # Filter date to create presence absence matrices for core and accessory genes
    get_pop_acc_pres_abs(input_dir, output_dir, acc_upper, acc_lower)
    get_pop_core_pres_abs(input_dir, output_dir, core_lower)

    print("Data filtered, presence absence matrices created for subpopulation of samples.")

    ## Get linkage matrices for CLARC analysis
    get_linkage_matrices(output_dir)

    print("Linkage matrices generated for the subpopulation accessory genes.")

    ## Perform blastn all v all comparison for accessory gene rep sequences
    subprocess.run(['bash', 'acccog_blastn.sh', output_dir])

    print("All vs. all nucleotide BLAST performed for the subpopulation accessory genes.")

    ## Perform eggnog functional annotation
    get_functional_groups(output_dir)

    print("Eggnog functional annotation complete. Only a little bit to go!")

    ## Perform CLARC cleaning of COGs
    clarc_cleaning(output_dir)

    print("CLARC finished re-defining COGs. Have fun with the results!")


if __name__ == "__main__":
    main()

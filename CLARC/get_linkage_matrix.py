#!/usr/bin/env python
# coding: utf-8

# # Get linkage matrices for the accessory genes in the specified subpopulation of samples

# This script calculates the necessary linkage matrices to run CLARC given a presence absence matrix of accessory genes in a given subpopulation of a more general pangenome analysis.
#
# For each COG pair, this script first does a pairwise count of: the number of times they both appear (p11), the number of times neither appears (p00), and the number of times that one appears but not the other (p01 and p10). Then it calculates a linkage metric by taking the log of the product log[p00*p11)/(p01*p10)], this is also a pairwise metric.
#
# So, this script outputs 5 files:
#
# - p00, p11, p10 and p01 matrices with counts for each state
# - log co-occurrence matrix

### Import necessary packages

import pandas as pd
import numpy as np
import scipy as scipy
import itertools
import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt
import fastcluster
import os
from joblib import Parallel, delayed

### Define function that calculates linkage matrices

def get_linkage_matrices(out_path):

    pres_abs_path = out_path+'/population_accessory_presence_absence.csv'

    # Read presence absence matrix
    pan_acc = pd.read_csv(pres_abs_path, index_col=0)

    # Set accession column as indices (since it's the isolate identifier)
    pan_acc.set_index('Accession', inplace=True)

    # Get list of COG names
    cog_names = list(pan_acc.columns.values)

    # Get number of accessory COGs (this will be the size of the matrices)
    cog_num = len(cog_names)

    # Define a function to calculate co-occurrence
    def calculate_occurrences(cog1, cog2, cog1_series, cog2_series):

        occur = (cog1_series.astype(str) + cog2_series.astype(str)).astype('category')

        # Count occurrences using value_counts
        p00 = (occur == '00').sum()
        p11 = (occur == '11').sum()
        p10 = (occur == '10').sum()
        p01 = (occur == '01').sum()

        return p00, p11, p10, p01

    # Create a blank matrix with n x n dimensions
    linkage_p00_result = pd.DataFrame(np.nan, index=cog_names, columns=cog_names)
    linkage_p11_result = pd.DataFrame(np.nan, index=cog_names, columns=cog_names)
    linkage_p10_result = pd.DataFrame(np.nan, index=cog_names, columns=cog_names)
    linkage_p01_result = pd.DataFrame(np.nan, index=cog_names, columns=cog_names)

    # Parallelize
    results = Parallel(n_jobs=-1)(
        delayed(calculate_occurrences)(cog1, cog2, pan_acc[cog1], pan_acc[cog2])
        for i, cog1 in enumerate(cog_names)
        for j, cog2 in enumerate(cog_names[i+1:], start=i+1)  # Start j from i+1 to avoid redundant calculations
    )

    # Assign values to result matrices
    index = 0
    for i, cog1 in enumerate(cog_names):
        for j, cog2 in enumerate(cog_names[i+1:], start=i+1):
            linkage_p00_result.at[cog1, cog2], linkage_p11_result.at[cog1, cog2], \
                linkage_p10_result.at[cog1, cog2], linkage_p01_result.at[cog1, cog2] = results[index]
            # Since matrices are symmetric, we can assign values to their mirrored positions
            linkage_p00_result.at[cog2, cog1], linkage_p11_result.at[cog2, cog1], \
                linkage_p10_result.at[cog2, cog1], linkage_p01_result.at[cog2, cog1] = results[index]
            index += 1

    # Convert dataframes to numpy arrays
    p11_array = linkage_p11_result.to_numpy()
    p00_array = linkage_p00_result.to_numpy()
    p01_array = linkage_p01_result.to_numpy()
    p10_array = linkage_p10_result.to_numpy()

    # Initialize the result matrix with NaN values
    linkage_logcoccur_result = np.full((cog_num, cog_num), np.nan)

    # Calculate log co-occurrence using vectorized operations
    p10p01 = p10_array * p01_array
    zero_indices = p10p01 == 0
    logl = np.log(np.divide(p11_array * p00_array, p10p01, where=~zero_indices))
    logl[zero_indices] = np.nan

    # Assign values to the result matrix
    linkage_logcoccur_result = pd.DataFrame(logl, index=cog_names, columns=cog_names)
    np.fill_diagonal(linkage_logcoccur_result.values, np.nan)

    # Create linkage folder to put the results there
    # Specify the directory path you want to create
    link_path = out_path+"/linkage"

    # Check if the directory already exists
    if not os.path.exists(link_path):
        # Create the directory
        os.makedirs(link_path)

    # Export results
    linkage_p00_result.to_csv(link_path+"/acc_p00_matrix.csv", index=True)
    linkage_p11_result.to_csv(link_path+"/acc_p11_matrix.csv", index=True)
    linkage_p10_result.to_csv(link_path+"/acc_p10_matrix.csv", index=True)
    linkage_p01_result.to_csv(link_path+"/acc_p01_matrix.csv", index=True)
    linkage_logcoccur_result.to_csv(link_path+"/acc_linkage_co-occur.csv", index=True)

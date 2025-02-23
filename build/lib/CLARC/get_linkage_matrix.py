#!/usr/bin/env python
# coding: utf-8

# # Get linkage matrices for the accessory genes in the specified subpopulation of samples

# This script calculates the necessary linkage matrices to run CLARC given a presence absence matrix of accessory genes in a given subpopulation of a more general pangenome analysis.
#
# For each COG pair, this script first does a pairwise count of: the number of times they both appear (p11), the number of times neither appears (p00), and the number of times that one appears but not the other (p01 and p10). Then it calculates a linkage metric by taking the log of the product log[p00*p11)/(p01*p10)], this is also a pairwise metric.
#
# So, this script outputs 4 files:
#
# - p00, p11, p10 and p01 matrices with counts for each state

### Import necessary packages

from joblib import Parallel, delayed
import pandas as pd
import numpy as np
import itertools
import os

# Starting from CLARC v1.1.3
#Revision list:

#- We don't need the linkage co-occur matrix. I've eliminated that part of the code.
#- Optimizing in a few ways:
      #- Instead of a loop, convert presence absence matrix to string and count co-occurance values
      #- Use numpy for this which is more efficient
      #- Because the presence absence only has 0s and 1s, we can use uint8 instead of int64 so save computing power

### Define function that calculates linkage matrices

def get_linkage_matrices(out_path, max_cores):

    max_cores = max_cores // 2

    def compute_pairwise_counts(i, j, pan_acc_str):
        """Instead of looping like we had previously been doing, we will compute the pairwise counts with numpy. This should be more efficient and improve CLARC runtime"""
        pairs = np.char.add(pan_acc_str[:, i], pan_acc_str[:, j])
        unique, counts = np.unique(pairs, return_counts=True)
        return i, j, dict(zip(unique, counts))

    pres_abs_path = os.path.join(out_path, 'population_accessory_presence_absence.csv')

    # Read presence absence matrix efficiently and convert to uint8 to save memory
    pan_acc = pd.read_csv(pres_abs_path, index_col=0).astype(np.uint8)

    # Convert dataframe to a NumPy string array for efficiency
    pan_acc_str = pan_acc.to_numpy(dtype='<U1')

    # Get column names
    columns = pan_acc.columns
    matrix_shape = (len(columns), len(columns))

    # Initialize matrices as NumPy arrays for faster computation
    P11_matrix = np.zeros(matrix_shape, dtype=int)
    P10_matrix = np.zeros(matrix_shape, dtype=int)
    P01_matrix = np.zeros(matrix_shape, dtype=int)
    P00_matrix = np.zeros(matrix_shape, dtype=int)

    # Convert column names to a list for indexing
    col_list = list(columns)
    col_index = {col: i for i, col in enumerate(col_list)}

    # Compute pairwise counts in parallel
    results = Parallel(n_jobs=max_cores, backend='loky')(delayed(compute_pairwise_counts)(i, j, pan_acc_str)
                                                   for i, j in itertools.combinations(range(len(columns)), 2))

    # Populate matrices with computed counts
    for i, j, count_dict in results:
        P11_matrix[i, j] = P11_matrix[j, i] = count_dict.get('11', 0)
        P10_matrix[i, j] = P10_matrix[j, i] = count_dict.get('10', 0)
        P01_matrix[i, j] = P01_matrix[j, i] = count_dict.get('01', 0)
        P00_matrix[i, j] = P00_matrix[j, i] = count_dict.get('00', 0)

    # Convert NumPy matrices to Pandas DataFrames for exporting
    P11_df = pd.DataFrame(P11_matrix, index=columns, columns=columns)
    P10_df = pd.DataFrame(P10_matrix, index=columns, columns=columns)
    P01_df = pd.DataFrame(P01_matrix, index=columns, columns=columns)
    P00_df = pd.DataFrame(P00_matrix, index=columns, columns=columns)

    # Create linkage folder
    link_path = os.path.join(out_path, "linkage")
    os.makedirs(link_path, exist_ok=True)

    # Export results
    P00_df.to_csv(os.path.join(link_path, "acc_p00_matrix.csv"), index=True)
    P11_df.to_csv(os.path.join(link_path, "acc_p11_matrix.csv"), index=True)
    P10_df.to_csv(os.path.join(link_path, "acc_p10_matrix.csv"), index=True)
    P01_df.to_csv(os.path.join(link_path, "acc_p01_matrix.csv"), index=True)

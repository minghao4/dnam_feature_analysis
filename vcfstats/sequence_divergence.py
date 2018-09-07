#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculating sequence divergence (Pi)
"""

import numpy as np
import pandas as pd
import sys

# Files
"""
Scaffold sizes file have been prefiltered for scaffolds above a certain length,
in our case, scaffolds above 100kb. A simple awk one-liner can be used for this
purpose:

(awk '{if ($2 > 100000) {print $0}}' canSat3_scaffold_sizes.tsv |
    sort -k1,1V) > sorted_filtered_scaffold_sizes.tsv

Usage:
    python filter_reformat_snps_pi_calc.py
        bin_size
        sorted_snps.tsv
        sorted_filtered_scaff_sizes.tsv
        output_folderk
"""

bin_size = int(sys.argv[1])
all_snps = sys.argv[2]
filtered_scaffold_sizes = sys.argv[3]
output_folder = sys.argv[4]


### PART ONE
# Identify SNPs on scaffolds greater than 100kb in length, discard others.
filtered_snps = pd.DataFrame(index = [0], columns = ["#Scaffold", "Position", "Allele_Frequency"])
snps_df = pd.read_csv(all_snps, sep = '\t')
freq_idx = snps_df.shape[1] - 1 # Index of final column
filtered_scaff_sizes = pd.read_csv(filtered_scaffold_sizes, sep = '\t')

currScaffold = ""
i = 0 # filtered_snps row index
status = False # Write status
for x in range(0, snps_df.shape[0]):
    row = snps_df.loc[x]
    scaffold_pos = row[0].split("_")
    scaffold = scaffold_pos[0]
    pos = scaffold_pos[1]
    freq = row[freq_idx]

    if scaffold != currScaffold:
        currScaffold = scaffold

        if scaffold in filtered_scaff_sizes.iloc[:,0].values:
            filtered_snps.loc[i] = [scaffold, pos, freq]
            i += 1
            status = True

        else:
            status = False

    elif status:
        filtered_snps.loc[i] = [scaffold, pos, freq]
        i += 1


### PART 2
# Calculate sequence divergence for each bin across EACH scaffold.

## Functions:

# Initialize output matrix
def init_matrix(num_bin, bin_size, header):
    # Pandas DataFrame of Numpy array
    # Sets midpoint of bin as label (e.g. bin 0-10000, labeled as 5000)
    curr_output = ((np.arange(num_bin * 2).reshape(num_bin, 2) + 1) * (bin_size / 2)).astype(float)
    curr_output[:, 1] *= 0 # Sets all values in second column (Pi) as 0.
    curr_output = pd.DataFrame(curr_output, columns = header)

    return curr_output

# Final pi calculation step.
def final_pi_calc(curr_output, bin_size):
    curr_output.iloc[:, 1] = curr_output.iloc[:, 1] * 2 / bin_size
    return curr_output

# Write output.
def write_output(scaff_name, curr_output):
    output_file = output_folder + "/" + scaff_name + "_pi.tsv"
    curr_output.to_csv(output_file, sep = '\t', index = False)

# Final pi calc and write output.
def final_calc_write(curr_output, bin_size, scaff_name):
    write_output(scaff_name, final_pi_calc(curr_output, bin_size))


## File I/O
header = ["#Distance", "Pi"]

# Current read scaffold values.
snp_idx = 0 # So we don't have to restart reading the filtered SNPs dataframe
for x in range(0, filtered_scaff_sizes.shape[0]):
    scaffold = filtered_scaff_sizes.loc[x]
    scaff_name = scaffold[0]
    num_bin = int(int(scaffold[1]) / bin_size) # Number of FULL bins, final partial bin excluded.
    curr_output = init_matrix(num_bin, bin_size, header)

    # Loop through SNP file.
    for j in range(snp_idx, filtered_snps.shape[0]):
        snp = filtered_snps.loc[j]
        scaff = snp[0] # Scaffold of SNP
        pos = int(snp[1]) # Position of SNP
        freq = float(snp[2]) # Allele frequency of SNP

        # If variant is on current scaffold
        if scaff == scaff_name:
            snp_idx += 1
            bin_idx = int(pos / bin_size)
            if (bin_idx < num_bin):
                curr_output.iloc[bin_idx, 1] += float((freq * (1 - freq) * 24 / 23))

        else:
            final_calc_write(curr_output, bin_size, scaff_name)
            break

final_calc_write(curr_output, bin_size, scaff_name)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
For parsing vcfs
"""

from pandas import DataFrame as df
import numpy as np
import os

# General variables
base_dict = {
    'A' : 1,
    'C' : 2,
    'G' : 3,
    'T' : 4,}

# General functions
def gt_status(gt):
    status = 0
    if gt == "0/1":
        status = 1

    return status

def allele_value(allele):
    val = 0
    for i in range(0, allele.length):
        base = allele[i]
        modifier = i + 1

        val += base_dict.get(base) * modifier

    return val

def variant_value(ref_allele, alt_allele):
    var_val = allele_value(alt_allele)
    if ref_allele.length != alt_allele.length:
        var_val -= allele_value(ref_allele)

    return var_val

def variant_dict_value(gt_status, ref_allele, alt_allele):
    alt_allele1_val = 0
    alt_allele2_val = 0

    if gt_status == 0:
        alt_allele1_val = variant_value(ref_allele, alt_allele)

    alt_allele2_val = variant_value(ref_allele, alt_allele)
    var_dict_entry = [alt_allele1_val, alt_allele2_val]

    return var_dict_entry


class VcfParser:
    def __init__(self, num_cultiv):
        self.variant_dict = {}
        self.curr_file = None
        self.num_cultiv = num_cultiv

    def read_curr_file(self, vcf_file_path):
        self.curr_file = df.read_csv(vcf_file_path, sep = '\t', index = False)
        self.curr_file.drop(columns = ["ID", "QUAL", "FILTER", "INFO", "FORMAT"])

    def process_curr_file(self, curr_cultiv_idx):
        for i in range(0, self.curr_file.shape[0]):
            row = self.curr_file.iloc[i]
            scaff_pos = row[0] + "_" + row[1]

            # Genotype, reference allele, and alternate allele, respectively.
            new_entry = variant_dict_value(row[4], row[2], row[3])

            if self.variant_dict.__contains__(scaff_pos):
                self.variant_dict.get(scaff_pos)[curr_cultiv_idx] = new_entry

            else:
                var_dict_entry = [[self.num_cultiv][2]]
                var_dict_entry[curr_cultiv_idx] = new_entry

                self.variant_dict.__setitem__(scaff_pos, var_dict_entry)


    def read_all_files(self, vcf_folder_path):
        ##TODO
        return



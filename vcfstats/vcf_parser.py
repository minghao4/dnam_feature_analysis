#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A parser for VCF files.
"""

from pandas import DataFrame as df
from pandas import read_csv
import numpy as np
import os

# General functions
def gt_status(gt):
    """
    Returns the genotype status (int).
    * 0 = Homozygous
    * 1 = Heterozygous
    """

    status = 0
    if gt == "0/1":
        status = 1

    return status

def variant_dict_value(gt_status, alt_allele, alt_freq):
    """
    Returns the variant dictionary entry (array) given the genotype status
    (int), reference allele (string), and alternate allele (string).
    """

    cultiv_allele1 = alt_allele
    cultiv_allele2 = alt_allele

    if gt_status == 0:
        cultiv_allele1 = "-"

    var_dict_entry = [cultiv_allele1, cultiv_allele2, alt_freq * 2]
    return var_dict_entry

def matrix_init(num_cultiv):
    blank = "-"
    m = [blank] * num_cultiv

    for i in range(num_cultiv):
        m[i] = [blank, blank, 0]

    return m

def set_header(cultiv_names):
    header = ["#Scaffold_Position"]
    for cultivar in cultiv_names:
        header.append(cultivar + "_A1")
        header.append(cultivar + "_A2")

    header.append("Allele_Frequency")
    return header

class VcfParser:
    def __init__(self):
        self.variant_dict = {}
        self.context_dict = {}
        self.curr_file = None
        self.cultiv_names = []
        self.var_out = None
        self.cxt_out = None

    def __read_curr_file(self, vcf_file_path):
        self.curr_file = read_csv(vcf_file_path, sep = '\t', index_col = None)
        self.curr_file = self.curr_file.drop(columns = ["N5-CONTEXT", "FORMAT"])

    def __process_curr_file(self, curr_cultiv_idx):
        for i in range(0, self.curr_file.shape[0]):
            row = self.curr_file.iloc[i]
            scaff_pos = row[0] + "_" + str(row[1])
            gt_sp_af = row[5].split(":")

            # Genotype, reference allele, alternate allele, and alt_freq,
            # respectively.
            new_entry = variant_dict_value(gt_status(gt_sp_af[0]), row[3], gt_sp_af[2])

            if self.variant_dict.__contains__(scaff_pos):
                self.variant_dict.get(scaff_pos)[curr_cultiv_idx] = new_entry

            else:
                var_dict_entry = matrix_init(len(self.cultiv_names))
                var_dict_entry[curr_cultiv_idx] = new_entry

                self.variant_dict.__setitem__(scaff_pos, var_dict_entry)
                self.context_dict.__setitem__(scaff_pos, row[4])

    def __set_var_out(self):
        header = set_header(self.cultiv_names)
        self.var_out = df(columns = header)
        i = 0

        for key, entries in self.variant_dict.items():
            data_row = [key]
            var_count = 0.0

            for entry in entries:
                data_row.append(entry[0])
                data_row.append(entry[1])
                var_count += float(entry[2])

            data_row.append(var_count)
            self.var_out.loc[i].append(var_count)
            i += 1

    def __set_cxt_out(self):
        header = ["#Scaffold_Position", "Context"]
        self.context_out = df(data = self.context_dict).T.reset_index()
        self.context_out.columns = header

    def __write_output(self, output_folder_path):
        print
        self.__set_var_out()
        self.__set_cxt_out()

        var_file_path = output_folder_path + "/master_variants.tsv"
        cxt_file_path = output_folder_path + "/variant_context.tsv"

        self.var_out.to_csv(var_file_path, sep = '\t', index = False)
        self.context_out.to_csv(cxt_file_path, sep = '\t', index = False)

    def parse_all_vcfs(self, vcf_folder_path):
        output_folder = False
        if vcf_folder_path.endswith("/"):
            vcf_folder_path = vcf_folder_path[:-1]

        for curr_cultiv_idx, filename in enumerate(sorted(os.listdir(vcf_folder_path))):
            vcf_file_path = vcf_folder_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                if filename == "output":
                    output_folder = True

                continue

            self.cultiv_names.append(filename.split("_")[0])
            self.__read_curr_file(vcf_file_path)
            self.__process_curr_file(curr_cultiv_idx)

        output_folder_path = vcf_folder_path + "/output"
        if not output_folder:
            os.mkdir(output_folder_path)

        self.__write_output(output_folder_path)

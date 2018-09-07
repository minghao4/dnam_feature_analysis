#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A parser for processed BISCUIT outputs.

This tool assumes the following format for the input file.
* A tab delimited file.
* Contains the following headers/columns:
    - Chromosome
    - Position
    - Reference Allele
    - Alternate Allele
    - 2 Nucleotide Context
    - 5 Nucleotide Context
    - Format (Genotype:Spread:Beta-Value/Alternate Allele Frequency)
    - Genotype:Spread:Beta-Value/Alternate Allele Frequency
"""

import numpy as np
import os
from pandas import DataFrame as df
from pandas import read_csv
from subprocess import run

# General functions
def gt_status(gt):
    """
    Defines and returns the genotype status.
    * 0 = Homozygous Ref
    * 1 = Heterozygous
    * 2 = Homozygous Variant

    Integer -> Integer
    """

    status = 0
    if gt == "0/1":
        status = 1

    elif gt == "1/1":
        status = 2

    return status


def variant_dict_value(gt_status, alt_allele, alt_freq):
    """
    Defines and returns the variant dictionary entry array for a single
    cultivar.

    Integer, String, Float -> List[Float, Float]
    """

    cultiv_allele_a = 0.0
    cultiv_allele_b = 0.0

    if alt_allele != ".":
        if gt_status < 2:
            cultiv_allele_b = float(alt_freq) * 2

        else:
            cultiv_allele_a = float(alt_freq)
            cultiv_allele_b = float(alt_freq)

    var_dict_entry = [cultiv_allele_a, cultiv_allele_b]
    return var_dict_entry


def matrix_init(num_cultiv):
    """
    Initializes and returns the variant dictionary entry array for all
    cultivars.

    Integer -> List[List[Float, Float]]
    """

    m = [0] * num_cultiv
    for i in range(num_cultiv):
        m[i] = [0, 0]

    return m


def set_var_header(cultiv_names):
    """
    Defines and returns the header for the variant output file.

    List[String] -> List[String]
    """

    header = ["#Scaffold_Position"]
    for cultivar in cultiv_names:
        header.append(cultivar + "_A1")
        header.append(cultivar + "_A2")

    header.append("Allele_Frequency")
    return header


def shell_sort(folder_path, file_name):
    """
    Runs the shell sort command (by version and ignores header) for a given
    file in a folder. Places sorted file in same folder.

    String, String -> None
    """

    old_filepath = str(folder_path + "/" + file_name)
    sorted_filepath = str(folder_path + "/" + "sorted_" + file_name)

    run(
        [str("(head -n 1 " + old_filepath + " && tail -n +2 " + old_filepath + ")")
         , "|", "sort", "-k1,1V", ">", sorted_filepath,

        ]

    )


class VcfParser:
    def __init__(self):
        self.variant_dict = {}
        self.context_dict = {}
        self.curr_file = None
        self.cultiv_names = []
        self.var_out = None
        self.cxt_out = None


    def __read_curr_file(self, vcf_file_path):
        """
        Reads the current file within the VcfParser object.

        VcfParser, String -> VcfParser
        """
        self.curr_file = read_csv(
                            vcf_file_path,
                            sep = '\t',
                            index_col = None,
                            engine = 'python',

                         )

        # TODO: Need to fix this line so it drops columns by index instead of by
        # name.
        self.curr_file = self.curr_file.drop(columns = ["5-CONTEXT", "FORMAT"])


    def __new_dict_entries(self, curr_cultiv_idx, new_entry, scaff_pos, row):
        """
        Initializes and sets a new dictionary entry for all variants
        """

        var_dict_entry = matrix_init(len(self.cultiv_names))
        var_dict_entry[curr_cultiv_idx] = new_entry

        self.variant_dict.__setitem__(scaff_pos, var_dict_entry)
        self.context_dict.__setitem__(scaff_pos, [row[4]])


    def __process_curr_file(self, curr_cultiv_idx):
        for i in range(self.curr_file.shape[0]):
            row = self.curr_file.iloc[i]
            scaff_pos = row[0] + "_" + str(row[1])
            gt_sp_af = row[5].split(":")

            new_entry = variant_dict_value(
                            gt_status(gt_sp_af[0]),
                            row[3],
                            gt_sp_af[2],

                        )

            if self.variant_dict.__contains__(scaff_pos):
                self.variant_dict.get(scaff_pos)[curr_cultiv_idx] = new_entry

            else:
                self.__new_dict_entries(
                    curr_cultiv_idx,
                    new_entry,
                    scaff_pos,
                    row,

                )

    def __set_var_out(self):
        header = set_var_header(self.cultiv_names)
        self.var_out = df(columns = header)
        i = 0

        for key, entries in self.variant_dict.items():
            data_row = []
            data_row.append(key)
            var_count = 0.0

            for entry in entries:
                data_row.append(entry[0])
                data_row.append(entry[1])
                var_count += float(entry[0])
                var_count += float(entry[1])

            data_row.append(float(var_count / len(self.cultiv_names)))
            self.var_out.loc[i] = data_row
            i += 1


    def __set_cxt_out(self):
        header = ["#Scaffold_Position", "Context"]
        self.context_out = df(data = self.context_dict).T.reset_index()
        self.context_out.columns = header


    def __write_output(self, output_folder_path):
        self.__set_var_out()
        self.__set_cxt_out()

        var_file_name = "master_variants.tsv"
        var_file_path = output_folder_path + "/" + var_file_name
        cxt_file_name = "variant_context.tsv"
        cxt_file_path = output_folder_path + "/" + cxt_file_name

        self.var_out.to_csv(var_file_path, sep = '\t', index = False)
        self.context_out.to_csv(cxt_file_path, sep = '\t', index = False)

        shell_sort(output_folder_path, var_file_name)
        shell_sort(output_folder_path, cxt_file_name)


    def parse_all_vcfs(self, vcf_folder_path):
        output_folder = False
        if vcf_folder_path.endswith("/"):
            vcf_folder_path = vcf_folder_path[:-1]

        folder = sorted(os.listdir(vcf_folder_path))

        for filename in folder:
            vcf_file_path = vcf_folder_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                if filename == "output":
                    output_folder = True

                continue

            self.cultiv_names.append(filename.split("_")[0])

        curr_cultiv_idx = 0
        for filename in folder:
            vcf_file_path = vcf_folder_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                continue

            self.__read_curr_file(vcf_file_path)
            self.__process_curr_file(curr_cultiv_idx)
            curr_cultiv_idx += 1

        output_folder_path = vcf_folder_path + "/output"
        if not output_folder:
            os.mkdir(output_folder_path)

        self.__write_output(output_folder_path)

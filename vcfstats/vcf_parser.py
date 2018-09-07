#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A parser for processed biscuit VCF outputs.

This tool assumes the following format for the input file.
* No metadata.
* A tab delimited file.
* Contains the following headers/columns in this order:
    - Chromosome
    - Position
    - Reference Allele
    - Alternate Allele
    - 2 Nucleotide Context
    - 5 Nucleotide Context
    - Format (Genotype:Spread:Beta-Value/Alternate Allele Frequency)
    - Genotype:Spread:Beta-Value/Alternate Allele Frequency

biscuit can be found here: https://github.com/zwdzwd/biscuit

Preprocessing may be implemented as part of this program later, but as of now
the biscuit VCF outputs require manual preprocessing prior to input into this
tool.
"""

import numpy as np
import os
import sys
from pandas import DataFrame as df
from pandas import read_csv
from subprocess import run, check_output
import subprocess


# General functions.
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


def var_dict_value(gt_status, alt_allele, alt_freq):
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


def set_var_matrix(num_cultiv):
    """
    Initializes and returns the variant dictionary entry array for all
    cultivars.

    Integer -> List[List[Float, Float]]
    """

    m = [0.0] * num_cultiv * 2
    # for i in range(num_cultiv):
    #     m[i] = [0, 0]

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


def shell_sort_sep(dir_path, file_name):
    """
    Runs the shell sort command (by version and ignores header) for a given
    file in a directory. Places sorted file in same directory.

    String, String -> None
    """

    old_filepath = str(dir_path + "/" + file_name)
    sorted_filepath = str(dir_path + "/" + "sorted_" + file_name)

    # Sorting.
    cmd = str("(head -n 1 " + old_filepath + \
              " && tail -n +2 " + old_filepath + \
              " | sort -k1,1V) > " + sorted_filepath)


    check_output(cmd, shell = True)

    # Removing the original and renaming the sorted file.
    try:
        cmd = str("rm " + old_filepath + " && mv " + sorted_filepath + " " + old_filepath)
        run(cmd, shell = True, check = True)

    except subprocess.CalledProcessError as exc:
        print("Status: FAIL\n", exc.returncode, "\n", exc.output)


    # Removing the sorted file.
    # check_output(str("rm " + old_filepath), shell = True)

    # # Separating #Scaffold_Position using awk and writing under original filename
    # try:
    #     paste_cmd = [
    #         "paste", "<(awk", "'{split($1,scaff_pos,\"_\");",
    #         "print scaff_pos[1]\"\t\"scaff_pos[2]}'",
    #         str(sorted_filepath + ")"),
    #         "<(cut", "-f", "2-",
    #         str(sorted_filepath + ")"),
    #         ">", old_filepath,

    #     ]

    #     # str("paste <(awk '{split($1, scaff_pos, \"_\"); \
    #     #                             print scaff_pos[1]\"\t\"scaff_pos[2]}' " + \
    #     #                         sorted_filepath + ") " + \
    #     #                     "<(cut -f 2- " + sorted_filepath + ") " + \
    #     #                     "> " + old_filepath

    #     #         )

    #     subprocess.run(paste_cmd, shell = True, check = True)

    # except subprocess.CalledProcessError as exc:
    #     print("Status: FAIL", exc.returncode, exc.output)
    #     sys.exit()



# Class specific functions.
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
        Reads the current file within the VcfParser object and sets it as a
        pandas DataFrame.

        VcfParser, String -> VcfParser
        """
        self.curr_file = read_csv(
                            vcf_file_path,
                            sep = '\t',
                            index_col = None,
                            engine = 'python',

                         )

        self.curr_file = self.curr_file.drop(
                            [self.curr_file.columns[5],
                                self.curr_file.columns[6]],
                            axis = 1,

                        )


    def __new_dict_entries(self, a1_idx, a2_idx, new_entry, scaff_pos, row):
        """
        Initializes and sets a new dictionary entry for all variants.

        VcfParser, Integer, Integer, List[Float, Float], String, List[8] ->
        VcfParser
        """

        var_dict_entry = set_var_matrix(len(self.cultiv_names))
        var_dict_entry[a1_idx] = new_entry[0]
        var_dict_entry[a2_idx] = new_entry[1]

        self.variant_dict.__setitem__(scaff_pos, var_dict_entry)
        self.context_dict.__setitem__(scaff_pos, [row[4], row[2], row[3]])


    def __process_curr_file(self, curr_cultiv_idx):
        """
        Processes the currently read file within the VcfParser object.

        VcfParser, Integer -> VcfParser
        """

        curr_cultiv_allele1_idx = curr_cultiv_idx * 2
        curr_cultiv_allele2_idx = curr_cultiv_allele1_idx + 1

        for i in range(self.curr_file.shape[0]):
            row = self.curr_file.iloc[i]
            scaff_pos = row[0] + "_" + str(row[1])
            gt_sp_af = row[5].split(":")

            new_entry = var_dict_value(
                            gt_status(gt_sp_af[0]),
                            row[3],
                            gt_sp_af[2],

                        )

            if self.variant_dict.__contains__(scaff_pos):
                self.variant_dict.get(scaff_pos)[curr_cultiv_allele1_idx] = \
                    new_entry[0]
                self.variant_dict.get(scaff_pos)[curr_cultiv_allele2_idx] = \
                    new_entry[1]

            else:
                self.__new_dict_entries(
                    curr_cultiv_allele1_idx,
                    curr_cultiv_allele2_idx,
                    new_entry,
                    scaff_pos,
                    row,

                )


    def __set_var_output(self):
        """
        Sets the variant output pandas DataFrame using data from the variant
        dictionary.

        VcfParser -> VcfParser
        """

        header = set_var_header(self.cultiv_names)
        self.var_out = df(data = self.variant_dict).T.reset_index()
        self.var_out["freq"] = 0.0
        self.var_out.columns = header
        self.var_out.loc[:]["Allele_Frequency"] = \
            self.var_out.sum(axis = 1) / len(self.cultiv_names)


    def __set_cxt_output(self):
        """
        Sets the context output pandas DataFrame using data from the context
        dictionary.

        VcfParser -> VcfParser
        """

        header = ["#Scaffold_Position", "Context", "Reference", "Alternate"]
        self.cxt_out = df(data = self.context_dict).T.reset_index()
        self.cxt_out.columns = header


    def __write_output(self, output_dir_path):
        """
        Writes the output files to disk.

        VcfParser, String -> VcfParser
        """

        self.__set_var_output()
        self.__set_cxt_output()

        var_file_name = "master_variants.tsv"
        var_file_path = output_dir_path + "/" + var_file_name
        cxt_file_name = "variant_context.tsv"
        cxt_file_path = output_dir_path + "/" + cxt_file_name

        self.var_out.to_csv(var_file_path, sep = '\t', index = False)
        self.cxt_out.to_csv(cxt_file_path, sep = '\t', index = False)

        shell_sort_sep(output_dir_path, var_file_name)
        shell_sort_sep(output_dir_path, cxt_file_name)


    def parse_all_vcfs(self, vcf_dir_path):
        """
        Main method. Parses all VCF files in a given directory and places final
        outputs in a dir/output directory. Output directory will be created if
        it doesn't already exist.

        VcfParser, String -> VcfParser.
        """

        output_dir_present = False
        if vcf_dir_path.endswith("/"):
            vcf_dir_path = vcf_dir_path[:-1]

        directory = sorted(os.listdir(vcf_dir_path))
        for filename in directory:
            vcf_file_path = vcf_dir_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                if filename == "output":
                    output_dir_present = True

                continue

            self.cultiv_names.append(filename.split("_")[0])

        curr_cultiv_idx = 0
        for filename in directory:
            vcf_file_path = vcf_dir_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                continue

            self.__read_curr_file(vcf_file_path)
            self.__process_curr_file(curr_cultiv_idx)
            curr_cultiv_idx += 1

        output_dir_path = vcf_dir_path + "/output"
        if not output_dir_present:
            os.mkdir(output_dir_path)

        self.__write_output(output_dir_path)

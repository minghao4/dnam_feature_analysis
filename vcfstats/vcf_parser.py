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

import os
from subprocess import run, check_output

from pandas import DataFrame as df
from pandas import read_csv


## General functions.

# Finds the methylation level and variation frequency for a VCF row within the
# Biscuit output VCF.
def meth_var_dictValues(alternate_allele, variant_frequency):
    """
    Defines and returns the methylation and variant dictionary entries for a
    single cultivar as an array.

    String, Float -> List[Float]
    """

    methylation_level = 1.0
    variation_frequency = 0.0
    if alternate_allele != ".":
        variation_frequency = variant_frequency

        if alternate_allele != "C" and alternate_allele != "G":
            methylation_level = 1.0 - variant_frequency

    return [methylation_level, variation_frequency]


# Sets the output file headers for the methylation level and variation frequency
# files.
#   - 0 = methylation file
#   - 1 = variation file
def set_out_dfHeaders(cultiv_names, output):
    """
    Defines and returns the header for the methylation and variation output
    files.

    List[String]m String -> List[String]
    """

    # TODO: make the header variable pass through from the UI
    header = ["#Scaffold_Position"]
    if output == "C":
        header += ["Context", "Reference", "Alternate"]

    else:
        header += cultiv_names
        if output == "V":
            header.append("Variation_Frequency")

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
    cmd = str("rm " + old_filepath + \
        " && mv " + sorted_filepath + " " + old_filepath)

    run(cmd, shell = True, check = True)


# Class specific functions.
class VcfParser:
    def __init__(self):
        self.methylation_dict = {}
        self.variation_dict = {}
        self.context_dict = {}
        self.curr_file = None
        self.cultiv_names = []
        self.meth_out = None
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


    def __add_new_dictEntries(
            self,
            newEntries_needed,
            curr_cultivIdx,
            new_dictEntries,
            scaff_pos,
            row,

        ):

        """
        Initializes and sets new methylation and variation dictionary entries
        for all VCF rows.

        VcfParser, String, Integer, List[Float], String, List[Object] ->
        VcfParser
        """

        if newEntries_needed == "B" or newEntries_needed == "M":
            meth_dictEntry = [0.0] * len(self.cultiv_names)
            meth_dictEntry[curr_cultivIdx] = new_dictEntries[0]
            self.methylation_dict.__setitem__(scaff_pos, meth_dictEntry)

        if newEntries_needed == "B" or newEntries_needed == "V":
            var_dictEntry = [0.0] * len(self.cultiv_names)
            var_dictEntry[curr_cultivIdx] = new_dictEntries[1]
            self.variation_dict.__setitem__(scaff_pos, var_dictEntry)

        self.context_dict.__setitem__(scaff_pos, [row[4], row[2], row[3]])


    def __process_curr_file(self, curr_cultivIdx):
        """
        Processes the currently read file within the VcfParser object.

        VcfParser, Integer -> VcfParser
        """

        for i in range(self.curr_file.shape[0]):
            row = self.curr_file.iloc[i]
            scaff_pos = row[0] + "_" + str(row[1])
            alternate_allele = row[3]
            variant_frequency = float(row[5].split(":")[2])

            new_dictEntries = meth_var_dictValues(
                                  alternate_allele,
                                  variant_frequency,

                              )

            newEntries_needed = "N"
            if self.methylation_dict.__contains__(scaff_pos):
                self.methylation_dict.get(scaff_pos)[curr_cultivIdx] = \
                    new_dictEntries[0]

            else:
                newEntries_needed = "M"

            if self.variation_dict.__contains__(scaff_pos):
                self.variation_dict.get(scaff_pos)[curr_cultivIdx] = \
                    new_dictEntries[1]

            else:
                if newEntries_needed == "M":
                    newEntries_needed = "B"

                else:
                    newEntries_needed = "V"

            if newEntries_needed != "N":
                self.__add_new_dictEntries(
                    newEntries_needed,
                    curr_cultivIdx,
                    new_dictEntries,
                    scaff_pos,
                    row,

                )


    def __set_meth_output(self):
        """

        VcfParser -> VcfParser
        """

        header = set_out_dfHeaders(self.cultiv_names, "M")
        self.meth_out = df(data = self.methylation_dict).T.reset_index()
        self.meth_out.columns = header


    def __set_var_output(self):
        """
        Sets the variation output pandas DataFrame using data from the variation
        dictionary.

        VcfParser -> VcfParser
        """

        header = set_out_dfHeaders(self.cultiv_names, "V")
        self.var_out = df(data = self.variation_dict).T.reset_index()
        self.var_out["freq"] = 0.0
        self.var_out.columns = header
        self.var_out.loc[:]["Variation_Frequency"] = \
            self.var_out.sum(axis = 1) / len(self.cultiv_names)


    def __set_cxt_output(self):
        """
        Sets the context output pandas DataFrame using data from the context
        dictionary.

        VcfParser -> VcfParser
        """

        header = set_out_dfHeaders(self.cultiv_names, "C")
        self.cxt_out = df(data = self.context_dict).T.reset_index()
        self.cxt_out.columns = header


    def __write_output(self, output_dir_path):
        """
        Writes the output files to disk.

        VcfParser, String -> VcfParser
        """

        self.__set_meth_output()
        self.__set_var_output()
        self.__set_cxt_output()


        meth_file_name = "methylation_levels.tsv"
        meth_file_path = output_dir_path + "/" + meth_file_name
        var_file_name = "variant_frequency.tsv"
        var_file_path = output_dir_path + "/" + var_file_name
        cxt_file_name = "variant_context.tsv"
        cxt_file_path = output_dir_path + "/" + cxt_file_name

        self.meth_out.to_csv(meth_file_path, sep = '\t', index = False)
        self.var_out.to_csv(var_file_path, sep = '\t', index = False)
        self.cxt_out.to_csv(cxt_file_path, sep = '\t', index = False)

        shell_sort_sep(output_dir_path, meth_file_name)
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

        curr_cultivIdx = 1
        for filename in directory:
            vcf_file_path = vcf_dir_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                continue

            self.__read_curr_file(vcf_file_path)
            self.__process_curr_file(curr_cultivIdx)

        output_dir_path = vcf_dir_path + "/output"
        if not output_dir_present:
            os.mkdir(output_dir_path)

        self.__write_output(output_dir_path)

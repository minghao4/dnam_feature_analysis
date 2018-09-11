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
# import sys
import timeit

from pandas import DataFrame as df
from pandas import read_csv


# TODO: move this to user interface.
start = timeit.default_timer()
# input_path = sys.argv[1]


## General functions.

# Finds the methylation level and variation frequency for a VCF row within the
# Biscuit output VCF.
def meth_var_dictValues(alternate_allele, variant_frequency):
    """
    Defines and returns the methylation and variation dictionary entries for a
    single cultivar as an array.

    String, Float -> List[Float]
    """

    methylation_level = 1.0
    variation_frequency = 0.0
    if alternate_allele != ".":
        variation_frequency = variant_frequency

        non_issue = ["C", "G", "N"]
        if alternate_allele not in non_issue:
            methylation_level = 1.0 - variant_frequency

    return [methylation_level, variation_frequency]


# Pipes a shell command out to sort outputs, remove original outputs, and
# rename sorted outputs.
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


    # Reads the VCF file given as a dataframe.
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

                            # Something is broken, not sure why I need this but
                            # it fixed things. Need to figure out virtualenv at
                            # some point.
                            engine = 'python',

                         )

        # Drops the 5-CONTEXT and FORMAT columns.
        # TODO: change parameters so that read_csv() doesn't read these columns
        # at all.
        self.curr_file = self.curr_file.drop(
                            [self.curr_file.columns[5],
                                self.curr_file.columns[6]],
                            axis = 1,

                        )


    # Adds a new dictionary entry if needed.
    #   - "B" = both methylation and variant dictionaries need new entries.
    #   - "M" = only methylation
    #   = "V" = only variation
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


    # For each variant in the VCF file currently being read, add to dictionaries.
    # Start a new key for each scaffold in the dictionaries.
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

            # This is probably not needed, wrote this in a rush. Methylation
            # and variation dictionaries should have the same scaffold positions
            # for each entry so for every instance both should need a new key
            # since they encounter a new scaffold at the same time.
            # TODO: remove this.
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


    # Sets the output file headers for the output files.
    def __set_out_dfHeaders(self, output):
        """
        Defines and returns the header for the output
        files.

        VcfParser, String -> VcfParser
        """

        # TODO: make the header variable pass through from the UI
        header = ["#Scaffold_Position"]
        if output == "C":
            header += ["Context", "Reference", "Alternate"]

        else:
            header += self.cultiv_names
            if output == "V":
                header.append("Variation_Frequency")

        return header


    def __set_meth_output(self):
        """
        Sets the methylation output pandas DataFrame using data from the
        variation dictionary.

        VcfParser -> VcfParser
        """

        header = self.__set_out_dfHeaders("M")
        self.meth_out = df(data = self.methylation_dict).T.reset_index()
        self.meth_out.columns = header


    def __set_var_output(self):
        """
        Sets the variation output pandas DataFrame using data from the variation
        dictionary.

        VcfParser -> VcfParser
        """

        header = self.__set_out_dfHeaders("V")
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

        header = self.__set_out_dfHeaders("C")
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


    # Main method incorporating everything.
    def parse_all_vcfs(self, vcf_dir_path):
        """
        Main method. Parses all VCF files in a given directory and places final
        outputs in a dir/output directory. Output directory will be created if
        it doesn't already exist.

        VcfParser, String -> VcfParser.
        """

        # TODO: move the stdout printing calls to a user interface file.
        print("\nStart.\n")
        output_dir_present = False
        if vcf_dir_path.endswith("/"):
            vcf_dir_path = vcf_dir_path[:-1]

        # Collecting the names of cultivars into a list.
        directory = sorted(os.listdir(vcf_dir_path))
        print("Creating list of cultivars:")
        for filename in directory:
            vcf_file_path = vcf_dir_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                if filename == "output":
                    output_dir_present = True

                continue

            cultivar = filename.split("_")[0]
            print("    - " + cultivar)
            self.cultiv_names.append(cultivar)

        # Reading and processing each of the files.
        print("\nProcessing all VCF files present...")
        curr_cultivIdx = 0
        for filename in directory:
            vcf_file_path = vcf_dir_path + "/" + filename
            if os.path.isdir(vcf_file_path):
                continue

            cultivar = filename.split("_")[0]
            print("\n" + cultivar + " VCF file:")
            print("Reading...")
            self.__read_curr_file(vcf_file_path)
            print("Processing...")
            self.__process_curr_file(curr_cultivIdx)
            curr_cultivIdx += 1

        # Output.
        print("Making output directory if not already present...\n")
        output_dir_path = vcf_dir_path + "/output"
        if not output_dir_present:
            os.mkdir(output_dir_path)

        print("\nWriting output files to " + output_dir_path)
        self.__write_output(output_dir_path)

        # Runtime.
        # TODO: move this to user interface as well.
        raw_runtime = timeit.default_timer() - start
        runtime = int(raw_runtime)
        hours = 0
        minutes = 0
        seconds = runtime % 60
        if runtime >= 3600:
            hours = int(runtime / 3600)
            minutes = int((runtime % 3600) / 60)

        elif runtime >= 60:
            minutes = int(runtime / 60)

        print("\n==========")
        print("VCF parsing complete.")
        print("Raw runtime: " + str(raw_runtime) + "s.")
        print("Program runtime: " + \
            str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s.")
        print("==========")


# parser = VcfParser()
# parser.parse_all_vcfs(input_path)

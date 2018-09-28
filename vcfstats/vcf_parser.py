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

from . import helpers

import os
from subprocess import run, check_output
# import sys
import timeit
from typing import List

import pandas as pd
from pandas import DataFrame as df


start_time = timeit.default_timer()
# vcf_file_path = sys.argv[1]


class VcfParser:
    def __init__(self):
        self.methylation_dict = {}
        self.variation_dict = {}
        self.current_file = None
        self.cultivar_names = []
        self.methylation_df = None
        self.variation_df = None


    # Reads the VCF file given as a dataframe.
    def __read_current_file(self, vcf_file_path: str) -> None:
        """
        Reads the current file within the VcfParser object and sets it as a
        pandas DataFrame.
        """
        #TODO: if broken, use engine = "python"
        self.current_file = pd.read_table(vcf_file_path)

        # Drops the CONTEXT, 5-CONTEXT and FORMAT columns.
        self.current_file = self.current_file.drop(
            ["CONTEXT", "5-CONTEXT", "FORMAT"], axis = 1
        )


    # Adds a new dictionary entry.
    def __add_new_dict_entries(
            self, cultivar_idx: int, scaffold_position: str, row: list,
            new_dict_entries: tuple(float, float)
        ) -> None:
        """
        Initializes and sets new methylation and variation dictionary entries
        for all VCF rows.
        """
        # Methylation dictionary entry.
        methylation_dict_entry = [0.0] * len(self.cultivar_names)
        methylation_dict_entry[cultivar_idx] = new_dict_entries[0]
        self.methylation_dict.__setitem__(
            scaffold_position, methylation_dict_entry
        )

        # Variation dictionary entry.
        variation_dict_entry = [0.0] * len(self.cultivar_names)
        variation_dict_entry[cultivar_idx] = new_dict_entries[1]
        self.variation_dict.__setitem__(scaffold_position, variation_dict_entry)


    # Finds the methylation level and variation frequency for a VCF row within the
    # Biscuit output VCF.
    @staticmethod
    def __row_dict_values(
            alternate_allele: str, variant_frequency: float
        ) -> tuple(float, float):
        """
        Defines and returns the methylation and variation dictionary entries
        for a single cultivar as an array.
        """
        # Initialize variables.
        methylation_level = 1.0
        variation_frequency = 0.0

        # Non-methylation variant.
        if alternate_allele != ".":
            variation_frequency = variant_frequency

            # Means methylated under BS-seq context.
            non_issue = ["C", "G", "N"]
            if alternate_allele not in non_issue:
                methylation_level = 1.0 - variant_frequency

        return (methylation_level, variation_frequency)


    # For each variant in the VCF file currently being read, add to dictionaries.
    # Start a new key for each scaffold in the dictionaries.
    def __process_current_file(self, cultivar_idx: int) -> None:
        """
        Processes the currently read file within the VcfParser object.
        """
        for i in range(self.current_file.shape[0]):
            row = self.current_file.iloc[i]
            scaffold_position = helpers.string_builder((
                row[0], "_", str(row[1])
            ))
            alternate_allele = row[3]
            variant_frequency = float(row[5].split(":")[2])
            new_dict_entries = self.__row_dict_values(
                alternate_allele, variant_frequency
            )

            new_entries_needed = False
            if self.methylation_dict.__contains__(scaffold_position):
                self.methylation_dict.get(scaffold_position)[cultivar_idx] = \
                    new_dict_entries[0]

            elif self.variation_dict.__contains__(scaffold_position):
                self.variation_dict.get(scaffold_position)[cultivar_idx] = \
                    new_dict_entries[1]

            else:
                new_entries_needed = True

            if new_entries_needed:
                self.__add_new_dict_entries(
                    cultivar_idx, scaffold_position, row, new_dict_entries
                )


    # Sets the output file headers for the output files.
    def __set_output_df_headers(self, output: str) -> List[str]:
        """
        Defines and returns the header for the output files.
        """
        header = ["#Scaffold_Position"] + self.cultivar_names
        if output == "variation":
            header += ["Variation_Frequency"]

        return header


    def __set_methylation_output_df(self) -> None:
        """
        Sets the methylation output pandas DataFrame using data from the
        methylation dictionary.
        """
        header = self.__set_output_df_headers("methylation")
        self.methylation_df = df(data = self.methylation_dict).T.reset_index()
        self.methylation_df.columns = header


    def __set_variation_output_df(self) -> None:
        """
        Sets the variation output pandas DataFrame using data from the variation
        dictionary.        """

        header = self.__set_output_df_headers("variation")
        self.variation_df = df(data = self.variation_dict).T.reset_index()
        self.variation_df["freq"] = 0.0
        self.variation_df.columns = header
        self.variation_df.loc[:]["Variation_Frequency"] = \
            self.variation_df.sum(axis = 1) / len(self.cultivar_names)


    # Pipes a shell command out to sort outputs, remove original outputs, and
    # rename sorted outputs.
    @staticmethod
    def __shell_sort_sep(dir_path: str, file_name: str) -> None:
        """
        Runs the shell sort command (by version and ignores header) for a given
        file in a directory. Places sorted file in same directory.
        """
        old_file_path = helpers.string_builder((dir_path, '/', file_name))
        sorted_file_path = helpers.string_builder((
            dir_path, '/', "sorted_", file_name
        ))

        # Sorting.
        cmd = helpers.string_builder((
            "(head -n 1 ", old_file_path, " && tail -n +2 ", old_file_path,
            " | sort -k1,1V) > ", sorted_file_path
        ))
        check_output(cmd, shell = True)

        # Removing the original and renaming the sorted file.
        cmd = helpers.string_builder((
            "rm ", old_file_path, " && mv ", sorted_file_path, ' ',
            old_file_path,
        ))
        run(cmd, shell = True, check = True)


    # Main method incorporating everything.
    def parse_all_vcfs(self, vcf_dir_path: str) -> None:
        """
        Main method. Parses all VCF files in a given directory and places final
        outputs in dir/output
        """
        # TODO: move the stdout printing calls to a user interface file.
        print()
        print("Start.")
        helpers.remove_trailing_slash([vcf_dir_path])

        # Collecting the names of cultivars into a list.
        directory = sorted(os.listdir(vcf_dir_path))
        print()
        print("Creating list of cultivars:")
        for file_name in directory:
            vcf_file_path = vcf_dir_path + "/" + file_name
            if os.path.isdir(vcf_file_path):
                continue

            cultivar = file_name.split("_")[0]
            print("    - " + cultivar)
            self.cultivar_names.append(cultivar)

        # Reading and processing each of the files.
        print()
        print("Processing all VCF files present...")
        current_cultivar_idx = 0
        for file_name in directory:
            vcf_file_path = vcf_dir_path + "/" + file_name
            if os.path.isdir(vcf_file_path):
                continue

            cultivar = file_name.split("_")[0]
            print()
            print(helpers.string_builder((cultivar, " VCF file:")))
            print("Reading...")
            self.__read_current_file(vcf_file_path)
            print("Processing...")
            self.__process_current_file(current_cultivar_idx)
            current_cultivar_idx += 1

        # Setting output dataframes.
        self.__set_methylation_output_df()
        self.__set_variation_output_df()

        # Writing to output
        output_dir_path = helpers.string_builder((vcf_dir_path, '/', "output"))
        methylation_file_name = "methylation_levels.tsv"
        helpers.write_output(
            self.methylation_df, methylation_file_name, output_dir_path
        )
        variation_file_name = "variant_frequency.tsv"
        helpers.write_output(
            self.variation_df, variation_file_name, output_dir_path
        )
        self.__shell_sort_sep(output_dir_path, methylation_file_name)
        self.__shell_sort_sep(output_dir_path, variation_file_name)

        # Runtime.
        # TODO: move to user interface.
        helpers.print_program_runtime("VCF parsing", start_time)


# parser = VcfParser()
# parser.parse_all_vcfs(vcf_file_path)

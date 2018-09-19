#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Methods for calculating DNAm sequence divergence (Pi)

The original formula for nucleotide diversity can be found here:
https://en.wikipedia.org/wiki/Nucleotide_diversity

This tool implements a slightly different formula. As we are working with
methylation sites, we are looking at the number of sites deviating from
normal methylation status per methylation site (within a pre-defined bin size).

The input parameters required are as follows:
* Bin size
* Sorted variation frequency TSV as per the output from the VCF Parser.
    - Scaffold_Position
    - Variation_Frequency

* Sorted scaffold sizes file.
* Output directory path.
"""

import os
# import sys
import timeit

import math
import numpy as np
from pandas import DataFrame as df
from pandas import read_csv


# TODO: move this to user interface.
start = timeit.default_timer()

# bin_width = int(sys.argv[1])
# var_file = sys.argv[2]
# scaff_sizes_file = sys.argv[3]
# output_dir_path = sys.argv[4]
# header_out = ["#Distance", "Pi"]


# Class specific methods.
class PiCalculator:
    def __init__(self):
        self.bin_size = 0
        self.var_df = None
        self.scaff_size_df = None
        self.curr_out_df = None


    # Reads the input variation file as a dataframe.
    def __set_in_df(self, bin_width, var_file, scaff_sizes_file):
        """
        Reads the variation file within the PiCalculator object and sets it as
        a pandas Dataframe.

        PiCalculator, Integer, String, String -> PiCalculator
        """

        self.bin_size = bin_width
        self.var_df = read_csv(
                        var_file,
                        sep = '\t',
                        usecols = ['#Scaffold_Position', 'Variation_Frequency']

                      )


        scaff_pos = self.var_df.columns[0]
        self.var_df[['Scaffold', 'Position']] = \
            self.var_df[scaff_pos].str.split('_', expand = True)

        # Reordering the columns.
        self.var_df = self.var_df.drop(scaff_pos, axis = 1)
        cols = self.var_df.columns.tolist()
        cols = cols[1:] + [cols[0]]
        self.var_df = self.var_df[cols]

        self.scaff_size_df = read_csv(scaff_sizes_file, sep = '\t')


    # Final part of the pi calculation.
    def __final_pi_calc(self, sites):
        """
        Performs the final pi calculation step.

        PiCalculator, Integer -> PiCalculator
        """
        divide = sites
        if sites == 0:
            divide = 1

        self.curr_out_df.iloc[:, 1] = self.curr_out_df.iloc[:, 1] * 2 / divide


    def __write_output(self, scaff_name, output_dir_path):
        """
        Writes the current output file to disk.

        PiCalculator, String, String -> PiCalculator
        """

        output_file = str(output_dir_path + "/" + scaff_name + "_pi_dnam.tsv")
        self.curr_out_df.to_csv(output_file, sep = '\t', index = False)


    def __final_calc_write(self, scaff_name, sites, output_dir_path):
        """
        Combined final pi calculation and output writing.

        PiCalculator, String, Integer, String -> PiCalculator
        """

        print("Final pi calculations...")
        self.__final_pi_calc(sites)
        print("Writing output.\n")
        self.__write_output(scaff_name, output_dir_path)


    # Loops through list of variants, calculates pi, and accumulates number of
    # sites. Writes the output for a scaffold upon hitting a variant on a
    # different scaffold. Returns a bookmark index and the final number of
    # methylation sites on this scaffold.
    def __read_var_df(
            self,
            var_df_bookmark,
            curr_scaff,
            num_bin,
            output_dir_path,

        ):

        """
        Processes the variation dataframe.

        PiCalculator, Integer, String, Integer, String -> PiCalculator,
        List[Integer, Boolean]
        """

        sites = 0
        for var_idx in range(var_df_bookmark, self.var_df.shape[0]):
            var = self.var_df.loc[var_idx]
            var_scaff = var[0]
            var_pos = int(var[1])
            var_freq = float(var[2])

            if var_scaff == curr_scaff:
                sites += 1
                var_df_bookmark += 1
                bin_idx = int(var_pos / self.bin_size)

                if bin_idx < num_bin:
                    # TODO: pipe in number of cultivars
                    pi_part = float(var_freq * (1 - var_freq) * 24 / 23)
                    self.curr_out_df.iloc[bin_idx, 1] += pi_part

            else:
                self.__final_calc_write(curr_scaff, sites, output_dir_path)
                break

        return [var_df_bookmark, sites]


    def __set_pi_matrix(self, num_bin, header, final_bin_lab):
        """
        Defines and returns the output pandas DataFrame for a scaffold.

        Integer, Integer, List[String], Float -> DataFrame
        """

        self.curr_out_df = ((np.arange(num_bin * 2).reshape(num_bin, 2) + 1) *
            (self.bin_size / 2)).astype(float)

        self.curr_out_df[:, 1] *= 0
        self.curr_out_df = df(self.curr_out_df, columns = header)

        out_dfRows = self.curr_out_df.shape[0]
        if out_dfRows > 1:
            self.curr_out_df.loc[out_dfRows - 1][0] = final_bin_lab


    # Loops through list of scaffolds, sets bins, calculates pi, and
    # accumulates number of sites. Writes the output for a scaffold upon
    # hitting a variant on a different scaffold. Returns a bookmark index and the
    # final number of methylation sites on this scaffold.
    def __read_scaff_df(self, header_out, output_dir_path):
        """
        Processes the scaffold list dataframe.

        PiCalculator, List[String] -> PiCalculator
        """

        var_df_bookmark = 0
        for scaff_idx in range(self.scaff_size_df.shape[0]):
            scaffold = self.scaff_size_df.iloc[scaff_idx]
            scaff_name = scaffold[0]
            print("Currently reading: " + scaff_name)

            scaff_size = int(scaffold[1])
            num_bin = math.ceil(scaff_size / self.bin_size)
            final_bin_lab = int((scaff_size % self.bin_size) / 2) + \
                (self.bin_size * (num_bin - 1))

            self.__set_pi_matrix(num_bin, header_out, final_bin_lab)
            varDfBookmark_sites = self.__read_var_df(
                                                  var_df_bookmark,
                                                  scaff_name,
                                                  num_bin,
                                                  output_dir_path,

                                              )

            var_df_bookmark= varDfBookmark_sites[0]
            sites = varDfBookmark_sites[1]

            if var_df_bookmark == self.var_df.shape[0]:
                self.__final_calc_write(scaff_name, sites, output_dir_path)

        # # Guard against 0
        # self.__final_calc_write(final_scaff, final_sites, output_dir_path)


    # Main method.
    def calculate_pi_all_scaffolds(
            self,
            bin_width,
            var_file,
            scaff_sizes_file,
            header_out,
            output_dir_path,

        ):

        """
        Main method.

        PiCalculator, Integer, String, String, List[String], String ->
        PiCalculator
        """

        print("\nStart.\n")
        paths = [var_file, scaff_sizes_file, output_dir_path]
        for path in paths:
            if path.endswith('/'):
                path = path[:-1]


        print("Creating output directory if it doesn't already exist...\n")
        if not os.path.isdir(output_dir_path):
            os.mkdir(output_dir_path)

        print("Setting variant dataframe...\n")
        self.__set_in_df(bin_width, var_file, scaff_sizes_file)
        print("Reading scaffold dataframe...\n")
        self.__read_scaff_df(header_out, output_dir_path)

        # Runtime.
        # TODO: move this to user interface as well.
        raw_runtime = timeit.default_timer() - start
        runtime = int(raw_runtime)
        hours = 0
        minutes = 0
        seconds = runtime % 60
        if runtime >= 3600:
            hours = runtime / 3600
            minutes = runtime % 3600 / 60

        elif runtime >= 60:
            minutes = runtime / 60

        print("\n==========")
        print("Sequence Divergence calculations complete.")
        print("Raw runtime: " + str(raw_runtime) + "s.")
        print("Program runtime: " + \
            str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s.")
        print("==========")


# pi_calc = PiCalculator()
# pi_calc.calculate_pi_all_scaffolds(
#     bin_width,
#     var_file,
#     scaff_sizes_file,
#     header_out,
#     output_dir_path,

# )

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
* Bin size [default 400].
* Sorted variants TSV as per the output from the VCF Parser.
    - Scaffold_Position
    - 2x Alleles for each cultivar
    - Allele_Frequency

* Sorted scaffold sizes file.
* Output directory path.
"""

import os
import sys

import math
import numpy as np
from pandas import DataFrame as df
from pandas import read_csv


# bin_width = int(sys.argv[1])
# var_file = sys.argv[2]
# scaff_sizes_file = sys.argv[3]
# output_dir_path = sys.argv[4]
# header_out = ["#Distance", "Pi"]


# General functions
def set_pi_matrix(num_bin, bin_size, header, final_bin_lab):
    """

    Integer, Integer, List[String], Float -> DataFrame
    """

    curr_out_df = ((np.arange(num_bin * 2).reshape(num_bin, 2) + 1) *
        (bin_size / 2)).astype(float)

    curr_out_df[:, 1] *= 0
    curr_out_df = df(curr_out_df, columns = header)
    if curr_out_df.shape[0] > 1:
        curr_out_df.loc[curr_out_df.shape[0] - 1][0] = final_bin_lab

    return curr_out_df


def final_pi_calc(curr_out_df, sites):
    """

    DataFrame, Integer -> DataFrame
    """
    print("before: ", curr_out_df.iloc[:, 1])
    curr_out_df.iloc[:, 1] = curr_out_df.iloc[:, 1] * 2 / sites
    print("after: ", curr_out_df.iloc[:, 1])
    return curr_out_df


# Class specific functions
class PiCalculator:
    def __init__(self):
        self.bin_size = 0
        self.var_df = df()
        self.scaff_size_df = df()
        self.curr_out_df = df()


    def __set_var_df(
            self,
            bin_width,
            var_file,
            scaff_sizes_file,
            output_dir_path,

        ):

        """

        PiCalculator, Integer, String, String, String -> PiCalculator
        """

        self.bin_size = bin_width
        self.var_df = read_csv(
                        var_file,
                        sep = '\t',
                        usecols = ['#Scaffold_Position', 'Allele_Frequency']

                      )


        self.var_df[['Scaffold', 'Position']] = \
            self.var_df['#Scaffold_Position'].str.split('_', expand = True)

        self.var_df = self.var_df.drop(['#Scaffold_Position'], axis = 1)
        cols = self.var_df.columns.tolist()
        cols = cols[1:] + [cols[0]]
        self.var_df = self.var_df[cols]
        self.scaff_size_df = read_csv(scaff_sizes_file, sep = '\t')


    def __write_output(self, scaff_name, output_dir_path):
        """

        PiCalculator, String, String -> PiCalculator
        """

        output_file = str(output_dir_path + "/" + scaff_name + "_pi_dnam.tsv")
        self.curr_out_df.to_csv(output_file, sep = '\t', index = False)


    def __final_calc_write(self, scaff_name, sites, output_dir_path):
        """

        PiCalculator, String, Integer, String -> PiCalculator
        """

        self.curr_out_df = final_pi_calc(self.curr_out_df, sites)
        self.__write_output(scaff_name, output_dir_path)


    def __read_var_df(
            self,
            var_df_bookmark,
            curr_scaff,
            num_bin,
            output_dir_path,

        ):

        """

        PiCalculator, Integer, String, Integer, String -> PiCalculator,
        List[Integer]
        """

        sites = 0
        for var_idx in range(var_df_bookmark, self.var_df.shape[0]):
            sites += 1
            var = self.var_df.loc[var_idx]
            var_scaff = var[0]
            var_pos = int(var[1])
            var_freq = float(var[2])
            write_status = False

            if var_scaff == curr_scaff:
                var_df_bookmark += 1
                bin_idx = int(var_pos / self.bin_size)

                if bin_idx < num_bin:
                    # TODO: pipe in number of cultivars
                    pi_part = float((var_freq * (1 - var_freq) * 24 / 23))

                    if pi_part > 0:
                        write_status = True

                    self.curr_out_df.iloc[bin_idx, 1] += pi_part

            else:
                if write_status:
                    print("Finishing touches and writing output...\n")
                    self.__final_calc_write(curr_scaff, sites, output_dir_path)

                break

        return [var_df_bookmark, sites]


    def __read_scaff_df(self, header_out, output_dir_path):
        """

        PiCalculator, List[String] -> PiCalculator
        """

        var_df_bookmark = 0
        final_scaff = ""
        final_sites = 0
        for scaff_idx in range(self.scaff_size_df.shape[0]):
            scaffold = self.scaff_size_df.iloc[scaff_idx]
            scaff_name = scaffold[0]
            scaff_size = int(scaffold[1])

            num_bin = math.ceil(scaff_size / self.bin_size)

            final_bin_lab = int((scaff_size % self.bin_size) / 2) + \
                (self.bin_size * (num_bin - 1))

            self.curr_out_df = set_pi_matrix(
                                   num_bin,
                                   self.bin_size,
                                   header_out,
                                   final_bin_lab,

                               )

            vdbm_s = self.__read_var_df(
                         var_df_bookmark,
                         scaff_name,
                         num_bin,
                         output_dir_path,

                     )

            var_df_bookmark= vdbm_s[0]
            final_scaff = scaff_name
            final_sites = vdbm_s[1]

        # Guard against 0
        self.__final_calc_write(final_scaff, final_sites + 1, output_dir_path)


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

        "Creating output directory if it doesn't already exist...\n"
        if not os.path.isdir(output_dir_path):
            os.mkdir(output_dir_path)

        "Setting variant dataframe...\n"
        self.__set_var_df(
            bin_width,
            var_file,
            scaff_sizes_file,
            output_dir_path,

        )

        "Reading scaffold dataframe...\n"
        self.__read_scaff_df(header_out, output_dir_path)


# pi_calc = PiCalculator()
# pi_calc.calulate_pi_all_scaffolds(
#     bin_width,
#     var_file,
#     scaff_sizes_file,
#     header_out,
#     output_dir_path,

# )

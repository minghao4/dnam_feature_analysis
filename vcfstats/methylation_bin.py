#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

import os
# import sys
import timeit

import numpy as np
from pandas import DataFrame as df
from pandas import read_csv
from scipy import stats as sps


start = timeit.default_timer()


def var_in_bin(curr_scaff, lb, ub, var_scaff, var_pos):
        """

        String, Integer, Integer, String, Integer -> Boolean
        """

        in_bin = False
        if var_scaff == curr_scaff:
            if var_pos >= lb and var_pos <= ub:
                in_bin = True

        return in_bin


class MethylationBinner:
    def __init__(self):
        self.var_methLvl_df = None
        self.bins_df = None


    def __set_in_df(self, bin_file, methLvl_file):
        """

        MethylationBinner, String, String -> MethylationBinner
        """

        self.bins_df = read_csv(bin_file, sep = "\t")
        self.var_methLvl_df = read_csv(methLvl_file, sep = "\t")

        scaff_pos = self.var_methLvl_df.columns[0]
        self.var_methLvl_df[['Scaffold', 'Position']] = \
            self.var_methLvl_df[scaff_pos].str.split('_', expand = True)

        # Reordering the columns.
        self.var_methLvl_df = self.var_methLvl_df.drop(scaff_pos, axis = 1)
        cols = self.var_methLvl_df.columns.tolist()
        cols = cols[(len(cols) - 2):] + cols[0:(len(cols) - 2)]
        self.var_methLvl_df = self.var_methLvl_df[cols]

        cultivs = cols[2:]
        for cultiv in cultivs:
            self.bins_df[cultiv] = 0


    def __read_var_methLvl_df(self, bin_idx, bookmark, curr_scaff, lb, ub):
        """

        MethylationBinner, Integer, String, Integer, Integer ->
        MethylationBinner, Integer
        """

        sites = 0
        for var_idx in range(bookmark, self.var_methLvl_df.shape[0]):
            var = self.var_methLvl_df.loc[var_idx]
            var_scaff = var[0]
            var_pos = int(var[1])

            if var_in_bin(curr_scaff, lb, ub, var_scaff, var_pos):
                print("Adding: " + var_scaff + " Position " + str(var_pos))
                sites += 1
                bookmark += 1
                self.bins_df.iloc[bin_idx, 2:] += var[2:]

            else:
                divide = 1
                if not sites == 0:
                    print("Averaging: " + curr_scaff)
                    divide = sites

                self.bins_df.iloc[bin_idx, 2:] /= divide
                break

        return bookmark



    def __read_bin_df(self):
        """

        MethylationBinner -> MethylationBinner
        """

        var_methLvl_df_bookmark = 0
        for bin_idx in range(self.bins_df.shape[0]):
            curr_bin = self.bins_df.loc[bin_idx]
            curr_scaff = curr_bin[0]
            curr_bin_lab = float(curr_bin[1])
            lb = 0
            ub = 0

            print()
            print("Reading: " + curr_scaff + " Bin " + str(curr_bin_lab))

            # TODO: fix this, 200bp is laziness because default bin is 400
            if curr_bin_lab % 200 == 0:
                lb = curr_bin_lab - 200 + 1
                ub = curr_bin_lab + 200

            else:
                lb = curr_bin_lab - curr_bin_lab % 200 + 1
                ub = curr_bin_lab + curr_bin_lab % 200

            var_methLvl_df_bookmark = self.__read_var_methLvl_df(
                                          bin_idx,
                                          var_methLvl_df_bookmark,
                                          curr_scaff,
                                          lb,
                                          ub,

                                      )

            if var_methLvl_df_bookmark == self.var_methLvl_df.shape[0]:
                break


    def __write_output(self, output_dir_path):
        """

        MethylationBinner, String -> MethylationBinner
        """

        output_file = str(output_dir_path + "/" + "methylation_bins.tsv")
        print("Creating output directory if it doesn't already exist...")
        print()
        if not os.path.isdir(output_dir_path):
            os.mkdir(output_dir_path)

        self.bins_df.to_csv(output_file, sep = '\t', index = False)


    def calculate_all_bin_methylation(
            self,
            bin_file,
            methLvl_file,
            output_dir_path,

        ):

        """
        Main method.

        MethylationBinner, String, String, String -> MethylationBinner
        """

        print()
        print("Start.")
        print()
        paths = [bin_file, methLvl_file, output_dir_path]
        for path in paths:
            if path.endswith('/'):
                path = path[:-1]

        print("Setting input dataframes...")
        print()
        self.__set_in_df(bin_file, methLvl_file)

        print("Calculating average methylation...")
        print()
        self.__read_bin_df()

        print("Writing to output...")
        print()
        self.__write_output(output_dir_path)

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

        print()
        print("==========")
        print("Methylation binning complete.")
        print("Raw runtime: " + str(raw_runtime) + "s.")
        print("Program runtime: " + \
            str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s.")
        print("==========")
        print()

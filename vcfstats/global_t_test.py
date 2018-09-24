#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from . import helpers

import os
# import sys
import timeit

from pandas import DataFrame as df
from pandas import read_csv
import scipy.stats as sps


start = timeit.default_timer()

class GlobalTTester:
    def __init__(self):
        self.lethbridge_df = None
        self.vegreville_df = None
        self.globalMeans_df = None
        self.out_df = None


    def __set_in_df(self, lethbridge_file, vegreville_file):
        """

        GlobalTTester, String, String -> GlobalTTester
        """

        self.lethbridge_df = read_csv(lethbridge_file, sep = '\t')
        self.vegreville_df = read_csv(vegreville_file, sep = '\t')

        mean_idx = self.lethbridge_df.columns[1:]
        mean_cols = ["Lethbridge", "Vegreville"]
        self.globalMeans_df = df(index = mean_idx, columns = mean_cols)

        out_idx = ["Global_T_Test"]
        out_cols = ["T_Statistic", "P-Value"]
        self.out_df = df(index = out_idx, columns = out_cols)

        self.globalMeans_df["Lethbridge"] = \
            self.lethbridge_df.mean(axis = 0).values

        self.globalMeans_df["Vegreville"] = \
            self.vegreville_df.mean(axis = 0).values


    def __global_t_test(self):
        """
        """

        model = sps.ttest_rel(
                    self.globalMeans_df["Vegreville"],
                    self.globalMeans_df["Lethbridge"],

                )

        self.out_df.iloc[0, 0:] = model[:]


    def global_t_testing(
            self,
            lethbridge_file,
            vegreville_file,
            output_dir_path,

        ):

        """
        Main Method.

        """

        self.__set_in_df(lethbridge_file, vegreville_file)
        self.__global_t_test()

        helpers.write_output(
            self.out_df,
            "Global_Methylation_TTest.tsv",
            output_dir_path,

        )

        helpers.print_program_runtime(
            "Global paired t-test calculations",
            start,

            )

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

import helpers

import os
# import sys
import timeit

from pandas import DataFrame as df
from pandas import read_csv
import scipy.stats as sps


start = timeit.default_timer()


class PairedTTester:
    def __init__(self):
        self.lethbridge_df = None
        self.vegreville_df = None
        self.bin_out_df = None
        self.cultiv_out_df = None



    def __set_in_df(self, lethbridge_file, vegreville_file):
        """

        PairedTTester, String, String -> PairedTTester
        """

        self.lethbridge_df = read_csv(lethbridge_file, sep = '\t')
        self.vegreville_df = read_csv(vegreville_file, sep = '\t')
        self.bin_out_df = self.lethbridge_df[self.lethbridge_df.columns[0:2]]
        self.bin_out_df['T-Statistic'] = 0
        self.bin_out_df['P-Value'] = 0
        self.cultiv_out_df = df(
                                 index = self.lethbridge_df.columns[2:],
                                 columns = self.bin_out_df.columns[2:],

                             )


    def __local_t_test(self):
        """
        """

        for bin_idx in range(self.bin_out_df.shape[0]):
            lethbridge_data = self.lethbridge_df.loc[bin_idx][2:]
            vegreville_data = self.vegreville_df.loc[bin_idx][2:]
            model = sps.ttest_rel(vegreville_data, lethbridge_data)

            self.bin_out_df.iloc[bin_idx, 2:] = model


    def __cultivar_t_test(self):
        """
        """

        for cultiv_idx in range(self.cultiv_out_df.shape[0]):
            cultiv = self.cultiv_out_df.iloc[:, cultiv_idx].name
            lethbridge_data = self.lethbridge_df[cultiv]
            vegreville_data = self.vegreville_df[cultiv]
            model = sps.ttest_rel(vegreville_data, lethbridge_data)

            self.cultiv_out_df.iloc[cultiv_idx] = model


    def paired_t_tests(self, lethbridge_file, vegreville_file, output_dir_path):
        """
        Main method.

        """

        self.__set_in_df(lethbridge_file, vegreville_file)
        self.__local_t_test()
        self.__cultivar_t_test()
        helpers.write_output(
            self.bin_out_df,
            "Cross_Variety_Methylation.tsv",
            output_dir_path,

        )

        helpers.write_output(
            self.bin_out_df,
            "Within_Variety_Methylation.tsv",
            output_dir_path,

        )

        helpers.print_program_runtime("Paired t-test calculations", start)

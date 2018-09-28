#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from . import helpers

import os
# import sys
import timeit

import pandas as pd
from pandas import DataFrame as df
import scipy.stats as sps


start_time = timeit.default_timer()
# lethbridge_file_path = sys.argv[1]
# vegreville_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]

class PairedTTester:
    def __init__(self):
        self.lethbridge_df = None
        self.vegreville_df = None
        self.bins_output_df = None
        self.cultivars_output_df = None


    def __set_dfs(
            self, lethbridge_file_path: str, vegreville_file_path: str
        ) -> None:
        """
        """
        self.lethbridge_df = pd.read_table(lethbridge_file_path)
        self.vegreville_df = pd.read_table(vegreville_file_path)

        self.bins_output_df = \
            self.lethbridge_df[self.lethbridge_df.columns[0:2]]
        self.bins_output_df['T-Statistic'] = 0
        self.bins_output_df['P-Value'] = 0

        self.cultivars_output_df = df(
            index = self.lethbridge_df.columns[2:],
            columns = self.bins_output_df.columns[2:]
        )


    def __local_t_test(self) -> None:
        """
        """
        for bin_idx in range(self.bins_output_df.shape[0]):
            lethbridge_data = self.lethbridge_df.loc[bin_idx][2:]
            vegreville_data = self.vegreville_df.loc[bin_idx][2:]
            model = sps.ttest_rel(vegreville_data, lethbridge_data)

            self.bins_output_df.iloc[bin_idx, 2:] = model[:]


    def __cultivar_t_test(self) -> None:
        """
        """
        for cultivar_idx in range(self.cultivars_output_df.shape[0]):
            cultivar = self.cultivars_output_df.iloc[:, cultivar_idx].name
            lethbridge_data = self.lethbridge_df[cultivar]
            vegreville_data = self.vegreville_df[cultivar]
            model = sps.ttest_rel(vegreville_data, lethbridge_data)

            self.cultivars_output_df.iloc[cultivar_idx] = model[:]


    def paired_t_tests(
            self, lethbridge_file_path: str, vegreville_file_path: str,
            output_dir_path: str
        ) -> None:
        """
        Main method.
        """
        self.__set_dfs(lethbridge_file_path, vegreville_file_path)
        self.__local_t_test()
        self.__cultivar_t_test()
        helpers.write_output(
            self.bins_output_df,
            "Cross_Variety_Methylation.tsv",
            output_dir_path,

        )

        helpers.write_output(
            self.cultivars_output_df,
            "Within_Variety_Methylation.tsv",
            output_dir_path,

        )

        helpers.print_program_runtime("Paired t-test calculations", start_time)


# ptt = PairedTTester()
# ptt.paired_t_tests(lethbridge_file_path, vegreville_file_path, output_dir_path)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import timeit, df, pd, sps
from . import helpers


start_time = timeit.default_timer()
# lethbridge_file_path = sys.argv[1]
# vegreville_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]


class GlobalTTester:
    def __init__(self):
        self.lethbridge_df = None
        self.vegreville_df = None
        self.global_means_df = None
        self.output_df = None


    def __set_dfs(
            self, lethbridge_file_path: str, vegreville_file_path: str
        ) -> None:
        """
        """
        self.lethbridge_df = pd.read_table(lethbridge_file_path)
        self.vegreville_df = pd.read_table(vegreville_file_path)
        self.lethbridge_df.drop(
            ["#Scaffold", "Bin_Label"], axis = 1, inplace = True
        )
        self.vegreville_df.drop(
            ["#Scaffold", "Bin_Label"], axis = 1, inplace = True
        )

        mean_idx = self.lethbridge_df.columns
        mean_cols = ["Lethbridge", "Vegreville"]
        self.global_means_df = df(index = mean_idx, columns = mean_cols)

        out_idx = ["Global_T_Test"]
        out_cols = ["T_Statistic", "P_Value", "Methylation_Ratio"]
        self.output_df = df(index = out_idx, columns = out_cols, data = 0)

        self.global_means_df["Lethbridge"] = \
            self.lethbridge_df.mean(axis = 0).values
        self.global_means_df["Vegreville"] = \
            self.vegreville_df.mean(axis = 0).values


    def __global_t_test(self) -> None:
        """
        """
        # Setting paired T-test results.
        model = sps.ttest_rel(
            self.global_means_df["Vegreville"],
            self.global_means_df["Lethbridge"]
        )
        self.output_df.iloc[0, 0:2] = model[:]

        # Setting methylation ratio.
        self.output_df.iloc[0, 2] = self.vegreville_df.sum().sum() / \
            self.lethbridge_df.sum().sum()
        self.output_df.index = ["global"]


    # Main method.
    def global_t_testing(
            self, lethbridge_file_path: str, vegreville_file_path: str,
            output_dir_path: str,
        ) -> None:
        """
        Main Method.
        """
        print("\nStart.")
        self.__set_dfs(lethbridge_file_path, vegreville_file_path)
        self.__global_t_test()

        helpers.write_output(
            self.output_df, "global_methylation_ttest.tsv", output_dir_path
        )

        helpers.print_program_runtime(
            "Global paired t-test calculations", start_time
        )


# gtt = GlobalTTester()
# gtt.global_t_testing(
#     lethbridge_file_path, vegreville_file_path, output_dir_path
# )

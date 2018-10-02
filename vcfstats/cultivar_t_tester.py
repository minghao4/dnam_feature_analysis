#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import timeit, math, df, pd, sps
from . import helpers


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
        self.bins_output_df["T-Statistic"] = 0
        self.bins_output_df["P-Value"] = 0
        self.bins_output_df["Significant?"] = False

        self.cultivars_output_df = df(
            data = 0,
            index = self.lethbridge_df.columns[2:],
            columns = self.bins_output_df.columns[2:]
        )
        self.cultivars_output_df["Significant?"] = False


    @staticmethod
    def __remove_nan(
            model: sps.stats.Ttest_relResult
        ) -> sps.stats.Ttest_relResult:
        """
        """
        model = \
            pd.Series(model).apply(lambda x: 0 if math.isnan(x) else x).values
        return model


    def __local_t_test(self) -> None:
        """
        """
        for bin_idx in range(self.bins_output_df.shape[0]):
            scaffold = self.bins_output_df.loc[bin_idx][0]
            bin_label = str(self.bins_output_df.loc[bin_idx][1])
            print(
                helpers.string_builder((
                    "T-Test: ", "Scaffold ", scaffold, " bin ", bin_label
                ))
            )
            lethbridge_data = self.lethbridge_df.loc[bin_idx][2:]
            vegreville_data = self.vegreville_df.loc[bin_idx][2:]
            model = self.__remove_nan(
                sps.ttest_rel(vegreville_data, lethbridge_data)
            )
            significant = helpers.significance(model.tolist())
            self.bins_output_df.iloc[bin_idx, 2:4] = model[:]
            self.bins_output_df.iloc[bin_idx, 4] = significant


    def __cultivar_t_test(self) -> None:
        """
        """
        for cultivar_idx in range(self.cultivars_output_df.shape[0]):
            cultivar = self.cultivars_output_df.iloc[cultivar_idx, :].name
            print(helpers.string_builder(("T-Test: ", cultivar)))
            lethbridge_data = self.lethbridge_df[cultivar]
            vegreville_data = self.vegreville_df[cultivar]
            # model = sps.ttest_rel(vegreville_data, lethbridge_data)
            model = self.__remove_nan(
                sps.ttest_rel(vegreville_data, lethbridge_data)
            )
            significant = helpers.significance(model.tolist())
            self.cultivars_output_df.iloc[cultivar_idx, 0:2] = model[:]
            self.cultivars_output_df.iloc[cultivar_idx, 2] = significant


    def paired_t_tests(
            self, lethbridge_file_path: str, vegreville_file_path: str,
            output_dir_path: str
        ) -> None:
        """
        Main method.
        """
        print()
        print("Start")
        helpers.remove_trailing_slash((
            lethbridge_file_path, vegreville_file_path, output_dir_path
        ))

        print()
        print("Setting input dataframes...")
        self.__set_dfs(lethbridge_file_path, vegreville_file_path)

        print()
        print("Performing local bin t-tests...")
        self.__local_t_test()

        print()
        print("Performing cultivar t-tests...")
        self.__cultivar_t_test()

        helpers.write_output(
            self.bins_output_df, "Cross_Variety_Methylation.tsv",
            output_dir_path
        )

        helpers.write_output(
            self.cultivars_output_df, "Within_Variety_Methylation.tsv",
            output_dir_path
        )

        helpers.print_program_runtime("Paired t-test calculations", start_time)


# ptt = PairedTTester()
# ptt.paired_t_tests(lethbridge_file_path, vegreville_file_path, output_dir_path)

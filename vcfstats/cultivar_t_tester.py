#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import List, timeit, math, df, pd, sps
from . import helpers

import multiprocessing
import sys


start_time = timeit.default_timer()
# lethbridge_file_path = sys.argv[1]
# vegreville_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]


class PairedTTesterInput:
    def __init__(self):
        self.lethbridge_df = None
        self.vegreville_df = None


    def set_input_dfs(
            self, lethbridge_file_path: str, vegreville_file_path: str
        ) -> None:
        """
        """
        self.lethbridge_df = pd.read_table(lethbridge_file_path)
        self.vegreville_df = pd.read_table(vegreville_file_path)


class LocalPairedTTestOutput:
    def __init__(self):
        self.bins_output_df = None


    def __set_output_df(self, methylation_input_df: pd.DataFrame) -> None:
        """
        """
        self.bins_output_df = methylation_input_df.iloc[:, 0:2]
        self.bins_output_df["T_Statistic"] = 0
        self.bins_output_df["P_Value"] = 1
        self.bins_output_df["Significant?"] = False


    def __local_t_test(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        for bin_idx in range(self.bins_output_df.shape[0]):
            bin_row = self.bins_output_df.loc[bin_idx]
            scaffold = bin_row[0]
            bin_label = str(bin_row[1])
            quick_search = helpers.string_builder((
                "T-Test: ", scaffold, '-', str(bin_label)
            ))

            lethbridge_data = \
                pd.to_numeric(lethbridge_input_df.loc[bin_idx][2:])
            vegreville_data = \
                pd.to_numeric(vegreville_input_df.loc[bin_idx][2:])
            if not lethbridge_data.equals(vegreville_data):
                model = sps.ttest_rel(vegreville_data, lethbridge_data)
                significant = helpers.significance(model)
                self.bins_output_df.iloc[bin_idx, 2:4] = model[:]
                self.bins_output_df.iloc[bin_idx, 4] = significant
                print(
                    helpers.string_builder((
                        '\n', "++++++++++", '\n', quick_search, '\n',
                        "T_Statistic: ", str(model[0]), '\n', "P_Value: ",
                        str(model[1]), '\n', '\n', "++++++++++", '\n'
                    ))
                )

            else:
                print(
                    helpers.string_builder((
                        '\n', scaffold, '-', str(bin_label),
                        " has been filtered...", '\n'
                    ))
                )


    def local_t_test_and_write(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame, output_dir_path: str
        ) -> None:
        """
        """
        print("Local t-test start.")
        sys.stdout = open("local_t_test_stdout.txt", 'w')
        print(
            helpers.string_builder((
                '\n', "-----**********-----", '\n', "Cross Variety T-Tests",
                '\n', "-----**********-----", '\n'
            ))
        )

        self.__set_output_df(lethbridge_input_df)
        self.__local_t_test(lethbridge_input_df, vegreville_input_df)

        helpers.write_output(
            self.bins_output_df, "cross_variety_methylation.tsv",
            output_dir_path
        )
        sys.stdout.close()


class CultivarPairedTTestOutput:
    def __init__(self):
        self.cultivars_output_df = None


    def __set_output_df(self, methylation_input_df: pd.DataFrame) -> None:
        """
        """
        self.cultivars_output_df = df(
            data = 0,
            index = methylation_input_df.columns[2:],
            columns = ["T_Statistic", "P_Value", "Significant?"]
        )
        self.cultivars_output_df["P_Value"] = 1
        self.cultivars_output_df["Significant?"] = False


    def __cultivar_t_test(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        """
        for cultivar_idx in range(self.cultivars_output_df.shape[0]):
            cultivar = self.cultivars_output_df.iloc[cultivar_idx, :].name
            print(
                helpers.string_builder((
                    '\n', "++++++++++", '\n', "T-Test: ", cultivar, '\n'
                ))
            )
            lethbridge_data = lethbridge_input_df[cultivar]
            vegreville_data = vegreville_input_df[cultivar]
            model = sps.ttest_rel(vegreville_data, lethbridge_data)
            # model = remove_inf(
            #     sps.ttest_rel(vegreville_data, lethbridge_data)
            # )
            significant = helpers.significance(model)
            self.cultivars_output_df.iloc[cultivar_idx, 0:2] = model[:]
            self.cultivars_output_df.iloc[cultivar_idx, 2] = significant
            print(
                helpers.string_builder((
                    "T_Statistic: ", str(model[0]), '\n', "P_Value: ",
                    str(model[1]), '\n', '\n', "++++++++++", '\n'
                ))
            )


    def cultivar_t_test_and_write(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame, output_dir_path: str
        ) -> None:
        """
        """
        print("Cultivar t-test start.")
        sys.stdout = open("cultivar_t_test_stdout.txt", 'w')
        print(
            helpers.string_builder((
                '\n', "-----**********-----", '\n', "Within Variety T-Tests",
                '\n', "-----**********-----", '\n'
            ))
        )

        self.__set_output_df(lethbridge_input_df)
        self.__cultivar_t_test(lethbridge_input_df, vegreville_input_df)

        helpers.write_output(
            self.cultivars_output_df, "within_variety_methylation.tsv",
            output_dir_path
        )
        sys.stdout.close()


def paired_t_tests(
        lethbridge_file_path: str, vegreville_file_path: str,
        output_dir_path: str
    ) -> None:
    """
    Main method.
    """
    print(helpers.string_builder(('\n', "Start.")))

    print(helpers.string_builder(('\n', "Setting dataframes...")))
    inputs = PairedTTesterInput()
    inputs.set_input_dfs(lethbridge_file_path, vegreville_file_path)

    print(
        helpers.string_builder((
            '\n', "Performing paired t-tests regression..."
        ))
    )
    local_output = LocalPairedTTestOutput()
    local_process = multiprocessing.Process(
        target = local_output.local_t_test_and_write,
        args = (
            inputs.lethbridge_df, inputs.vegreville_df, output_dir_path
        )
    )
    cultivar_output = CultivarPairedTTestOutput()
    cultivar_process = multiprocessing.Process(
        target = cultivar_output.cultivar_t_test_and_write,
        args = (
            inputs.lethbridge_df, inputs.vegreville_df, output_dir_path
        )
    )

    local_process.start()
    cultivar_process.start()
    local_process.join()
    cultivar_process.join()

    helpers.print_program_runtime("Paired t-test calculations", start_time)


# paired_t_tests(lethbridge_file_path, vegreville_file_path, output_dir_path)

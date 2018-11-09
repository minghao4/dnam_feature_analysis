#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import List, multiprocessing, sys, timeit, math, pd, df, sps
from . import helpers


class PairedTTesterInput:
    def __init__(self) -> None:
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
    def __init__(self) -> None:
        self.bins_output_df = None


    def __set_output_df(self, methylation_input_df: pd.DataFrame) -> None:
        """
        """
        self.bins_output_df = methylation_input_df.iloc[:, 0:2]

        # Default values.
        self.bins_output_df["T_Statistic"] = 0
        self.bins_output_df["P_Value"] = 1
        self.bins_output_df["Methylation_Ratio"] = 1
        self.bins_output_df["Significant?"] = False


    def __filter_and_test(
            self, lethbridge_data: str, vegreville_data: str, bin_idx: int,
            scaffold: str, bin_label: str, quick_search: str
        ) -> None:
        """
        """
        if not lethbridge_data.equals(vegreville_data):
            model = sps.ttest_rel(vegreville_data, lethbridge_data)
            significant = helpers.significance(model)
            methylation_ratio = vegreville_data.sum() / \
                lethbridge_data.sum()

            self.bins_output_df.iloc[bin_idx, 2:4] = model[:]
            self.bins_output_df.iloc[bin_idx, 4] = methylation_ratio
            self.bins_output_df.iloc[bin_idx, 5] = significant

            wrapping_flair = helpers.string_builder(('\n', '+' * 10, '\n'))
            print(helpers.string_builder((
                wrapping_flair, quick_search, "\nT_Statistic: ",
                str(model[0]), "\nP_Value: ", str(model[1]),
                "\nMethylation_Ratio: ", methylation_ratio,
                wrapping_flair
            )))

        else:
            print(helpers.string_builder((
                    '\n', scaffold, '-', str(bin_label),
                    " has been filtered...\n"
            )))


    def __iter_bins(
            self, bin_row: pd.Series, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        """
        bin_idx = bin_row.name
        scaffold = bin_row[0]
        bin_label = str(bin_row[1])
        quick_search = helpers.string_builder((
            "T-Test: ", scaffold, '-', str(bin_label)
        ))

        lethbridge_data = pd.to_numeric(
            lethbridge_input_df.loc[bin_idx][2:]
        )
        vegreville_data = pd.to_numeric(
            vegreville_input_df.loc[bin_idx][2:]
        )

        self.__filter_and_test(
            lethbridge_data, vegreville_data, bin_idx, scaffold, bin_label,
            quick_search
        )



    def __local_t_test(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        """
        tmp = self.bins_output_df.apply(
            self.__iter_bins, axis = 1,
            args = (lethbridge_input_df, vegreville_input_df)
        )
        del tmp
        # for bin_idx in range(self.bins_output_df.shape[0]):
        #     bin_row = self.bins_output_df.loc[bin_idx]
        #     scaffold = bin_row[0]
        #     bin_label = str(bin_row[1])
        #     quick_search = helpers.string_builder((
        #         "T-Test: ", scaffold, '-', str(bin_label)
        #     ))

        #     lethbridge_data = pd.to_numeric(
        #         lethbridge_input_df.loc[bin_idx][2:]
        #     )
        #     vegreville_data = pd.to_numeric(
        #         vegreville_input_df.loc[bin_idx][2:]
        #     )

        #     if not lethbridge_data.equals(vegreville_data):
        #         model = sps.ttest_rel(vegreville_data, lethbridge_data)
        #         significant = helpers.significance(model)
        #         methylation_ratio = vegreville_data.sum() / \
        #             lethbridge_data.sum()

        #         self.bins_output_df.iloc[bin_idx, 2:4] = model[:]
        #         self.bins_output_df.iloc[bin_idx, 4] = methylation_ratio
        #         self.bins_output_df.iloc[bin_idx, 5] = significant

        #         wrapping_flair = helpers.string_builder(('\n', '+' * 10, '\n'))
        #         print(helpers.string_builder((
        #             wrapping_flair, quick_search, "\nT_Statistic: ",
        #             str(model[0]), "\nP_Value: ", str(model[1]),
        #             "\nMethylation_Ratio: ", methylation_ratio,
        #             wrapping_flair
        #         )))

        #     else:
        #         print(helpers.string_builder((
        #             '\n', scaffold, '-', str(bin_label),
        #             " has been filtered...\n"
        #         )))


    def local_t_test_and_write(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame, output_dir_path: str
        ) -> None:
        """
        """
        print("Local t-test start.")

        # Prints stdout to separate file.
        sys.stdout = open("local_t_test_stdout.txt", 'w')
        wrapping_flair = helpers.string_builder((
            '\n', '-' * 5, '*' * 10, '-' * 5, '\n'
        ))
        print(helpers.string_builder((
            wrapping_flair, "Cross Variety T-Tests", wrapping_flair
        )))

        self.__set_output_df(lethbridge_input_df)
        self.__local_t_test(lethbridge_input_df, vegreville_input_df)
        helpers.write_output(
            output_df = self.bins_output_df,
            output_file_name = "cross_variety_methylation_ttest.tsv",
            output_dir_path = output_dir_path
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
            columns = [
                "Cultivar", "T_Statistic", "P_Value", "Methylation_Ratio",
                "Significant?"
            ]
        )

        # Defaults values.
        self.cultivars_output_df["Cultivar"] = methylation_input_df.columns[2:],
        self.cultivars_output_df["P_Value"] = 1
        self.cultivars_output_df["Methylation_Ratio"] = 1
        self.cultivars_output_df["Significant?"] = False


    def __filter_and_test(
            self, lethbridge_data: str, vegreville_data: str, cultivar_idx: int,
            cultivar: str
        ) -> None:
        """
        """
        model = sps.ttest_rel(vegreville_data, lethbridge_data)
        significant = helpers.significance(model)
        methylation_ratio = vegreville_data.sum() / \
            lethbridge_data.sum()

        self.cultivars_output_df.iloc[cultivar_idx, 1:3] = model[:]
        self.cultivars_output_df.iloc[cultivar_idx, 3] = methylation_ratio
        self.cultivars_output_df.iloc[cultivar_idx, 4] = significant

        wrapping_flair = helpers.string_builder(('\n', '+' * 10, '\n'))
        print(helpers.string_builder((
            "T_Statistic: ", str(model[0]), "\nP_Value: ", str(model[1]),
            "\nMethylation_Ratio: ", methylation_ratio, wrapping_flair
        )))


    def __iter_bins(
            self, cultivar_row: pd.Series, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        """
        cultivar_idx = cultivar_row.name
        cultivar = cultivar_row[0]
        print(helpers.string_builder((
            "\n++++++++++\n", "T-Test: ", cultivar, '\n'
        )))

        lethbridge_data = lethbridge_input_df[cultivar]
        vegreville_data = vegreville_input_df[cultivar]

        self.__filter_and_test(
            lethbridge_data, vegreville_data, cultivar_idx, cultivar
        )


    def __cultivar_t_test(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        """
        tmp = self.cultivars_output_df.apply(
            self.__iter_bins, axis = 1,
            args = (lethbridge_input_df, vegreville_input_df)
        )
        del tmp
        # for cultivar_idx in range(self.cultivars_output_df.shape[0]):
        #     cultivar = self.cultivars_output_df.iloc[cultivar_idx, :].name
        #     print(
        #         helpers.string_builder((
        #             '\n', "++++++++++", '\n', "T-Test: ", cultivar, '\n'
        #         ))
        #     )
        #     lethbridge_data = lethbridge_input_df[cultivar]
        #     vegreville_data = vegreville_input_df[cultivar]
        #     model = sps.ttest_rel(vegreville_data, lethbridge_data)
        #     significant = helpers.significance(model)
        #     self.cultivars_output_df.iloc[cultivar_idx, 0:2] = model[:]
        #     methylation_ratio = vegreville_data.sum() / \
        #             lethbridge_data.sum()
        #     self.cultivars_output_df.iloc[cultivar_idx, 2] = methylation_ratio
        #     self.cultivars_output_df.iloc[cultivar_idx, 3] = significant
        #     print(
        #         helpers.string_builder((
        #             "T_Statistic: ", str(model[0]), '\n', "P_Value: ",
        #             str(model[1]), '\n', "Methylation_Ratio: ",
        #             methylation_ratio, '\n', '\n', "++++++++++", '\n'
        #         ))
        #     )


    def cultivar_t_test_and_write(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame, output_dir_path: str
        ) -> None:
        """
        """
        print("Cultivar t-test start.")

        # Prints stdout to separate file.
        sys.stdout = open("cultivar_t_test_stdout.txt", 'w')
        wrapping_flair = helpers.string_builder((
            '\n', '-' * 5, '*' * 10, '-' * 5, '\n'
        ))
        print(helpers.string_builder((
            wrapping_flair, "Within Variety T-Tests", wrapping_flair
        )))

        self.__set_output_df(lethbridge_input_df)
        self.__cultivar_t_test(lethbridge_input_df, vegreville_input_df)

        helpers.write_output(
            output_df = self.cultivars_output_df,
            output_file_name = "within_variety_methylation_ttest.tsv",
            output_dir_path = output_dir_path
        )
        sys.stdout.close()


# Main method.
def paired_t_tests(
        lethbridge_file_path: str, vegreville_file_path: str,
        output_dir_path: str
    ) -> None:
    """
    Main method.
    """
    start_time = timeit.default_timer() # Initialize starting time.
    lethbridge_file_path, vegreville_file_path, output_dir_path = \
        helpers.remove_trailing_slash((
            lethbridge_file_path, vegreville_file_path, output_dir_path
        ))

    print("\nStart\nSetting dataframes...")  # Initialize input object.
    inputs = PairedTTesterInput()
    inputs.set_input_dfs(lethbridge_file_path, vegreville_file_path)

    print("\nPerforming paired t-tests regression...") # Initiate output objects.
    local_output = LocalPairedTTestOutput()
    cultivar_output = CultivarPairedTTestOutput()

    # Initialize processes for local and cultivar paired t-test objects.
    local_process = multiprocessing.Process(
        target = local_output.local_t_test_and_write,
        args = (
            inputs.lethbridge_df, inputs.vegreville_df, output_dir_path
        )
    )
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

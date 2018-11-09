#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Objective: perform cross-cultivar, within-cultivar, and global paired T-tests
between Vegreville and Lethbridge data.

Cross-cultivar: across all cultivar methylation data at a certain bin.
Within-cultivar: across all bin methylation data at a given cultivar.
Global: comparing global cultivar means across all cultivars.

Inputs:
- Lethbridge binned methylation TSV file path.
- Vegreville binned methylation TSV file path.
- Output directory path.

Outputs:
- Log TXT stdout files.
- TSV file holding paired T-test results (t-value, p-value, methylation ratio,
  nominal significance).

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
        Set input dataframes from given Lethbridge and Vegreville file paths.
        """
        self.lethbridge_df = pd.read_table(lethbridge_file_path)
        self.vegreville_df = pd.read_table(vegreville_file_path)


class LocalPairedTTestOutput:
    def __init__(self) -> None:
        self.bins_output_df = None


    def __set_output_df(self, methylation_input_df: pd.DataFrame) -> None:
        """
        Set output dataframes given an input dataframe's Scaffold and Position
        column data.
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
        Filter out bins where Lethbridge and Vegreville data is identical
        across all cultivars. If not identical, perform paired T-test.
        """
        # If the rows are equal, paired t-test cannot be performed.
        if not lethbridge_data.equals(vegreville_data):
            # Vegreville always listed first.
            model = sps.ttest_rel(vegreville_data, lethbridge_data)
            significant = helpers.significance(model) # bool
            methylation_ratio = vegreville_data.sum() / \
                (lethbridge_data.sum() + 0.01) # Avoid dividing by zero.

            self.bins_output_df.iloc[bin_idx, 2:4] = model[:]
            self.bins_output_df.iloc[bin_idx, 4] = methylation_ratio
            self.bins_output_df.iloc[bin_idx, 5] = significant

            # Printing to stdout to log findings.
            wrapping_flair = helpers.string_builder(('\n', '+' * 10, '\n'))
            print(helpers.string_builder((
                wrapping_flair, quick_search, "\nT_Statistic: ",
                str(model[0]), "\nP_Value: ", str(model[1]),
                "\nMethylation_Ratio: ", methylation_ratio,
                wrapping_flair
            )))

        else:
            # Printing to stdout to log findings.
            print(helpers.string_builder((
                    '\n', scaffold, '-', str(bin_label),
                    " has been filtered...\n"
            )))


    def __iter_bins(
            self, bin_row: pd.Series, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        Iterate through bins and perform paired T-tests between Vegreville and
        Lethbridge bin data if the bin has not been filtered.
        """
        bin_idx = bin_row.name
        scaffold = bin_row[0]
        bin_label = str(bin_row[1])

        # Quick search string in stdout log file.
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


    # def __local_t_test(
    #         self, lethbridge_input_df: pd.DataFrame,
    #         vegreville_input_df: pd.DataFrame
    #     ) -> None:
    #     """
    #     Perform cross-cultivar paired T-tests.
    #     """
    #     tmp = self.bins_output_df.apply(
    #         self.__iter_bins, axis = 1,
    #         args = (lethbridge_input_df, vegreville_input_df)
    #     )
    #     del tmp
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
        Perform cross-cultivar paired T-tests and save the output dataframe to
        a file.
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

        # Cross-cultivar paired T-tests.
        self.__set_output_df(lethbridge_input_df)
        tmp = self.bins_output_df.apply(
            self.__iter_bins, axis = 1,
            args = (lethbridge_input_df, vegreville_input_df)
        )
        del tmp
        helpers.write_output(
            output_df = self.bins_output_df,
            output_file_name = "cross_variety_methylation_ttest.tsv",
            output_dir_path = output_dir_path
        )
        sys.stdout.close()


class CultivarPairedTTestOutput:
    def __init__(self) -> None:
        self.cultivars_output_df = None


    def __set_output_df(self, methylation_input_df: pd.DataFrame) -> None:
        """
        Set output dataframes given an input dataframe's cultivars.
        """
        self.cultivars_output_df = df(
            data = 0,
            columns = [
                "Cultivar", "T_Statistic", "P_Value", "Methylation_Ratio",
                "Significant?"
            ]
        )

        # Default values.
        self.cultivars_output_df["Cultivar"] = methylation_input_df.columns[2:],
        self.cultivars_output_df["P_Value"] = 1
        self.cultivars_output_df["Methylation_Ratio"] = 1
        self.cultivars_output_df["Significant?"] = False


    def __cultivar_paired_t_test(
            self, lethbridge_data: str, vegreville_data: str, cultivar_idx: int,
            cultivar: str
        ) -> None:
        """
        Performed paired T-test between cultivar data.
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


    def __iter_cultivars(
            self, cultivar_row: pd.Series, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        Iterate through cultivars and perform paired T-tests between Vegreville
        and Lethbridge cultivar data.
        """
        cultivar_idx = cultivar_row.name
        cultivar = cultivar_row[0]
        print(helpers.string_builder((
            "\n++++++++++\n", "T-Test: ", cultivar, '\n'
        )))

        lethbridge_data = lethbridge_input_df[cultivar]
        vegreville_data = vegreville_input_df[cultivar]

        self.__cultivar_paired_t_test(
            lethbridge_data, vegreville_data, cultivar_idx, cultivar
        )


    # def __cultivar_t_test(
    #         self, lethbridge_input_df: pd.DataFrame,
    #         vegreville_input_df: pd.DataFrame
    #     ) -> None:
    #     """
    #     """
    #     tmp = self.cultivars_output_df.apply(
    #         self.__iter_cultivars, axis = 1,
    #         args = (lethbridge_input_df, vegreville_input_df)
    #     )
    #     del tmp
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
        Perform cross-cultivar paired T-tests and save output to a file.
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
        tmp = self.cultivars_output_df.apply(
            self.__iter_cultivars, axis = 1,
            args = (lethbridge_input_df, vegreville_input_df)
        )
        del tmp
        helpers.write_output(
            output_df = self.cultivars_output_df,
            output_file_name = "within_variety_methylation_ttest.tsv",
            output_dir_path = output_dir_path
        )
        sys.stdout.close()


class GlobalPairedTTestOutput:
    def __init__(self) -> None:
        self.global_means_df = None
        self.global_output_df = None


    def __set_output_dfs(self, cultivars: pd.Index) -> None:
        """
        """
        self.global_means_df = df(
            index = cultivars, columns = ["Lethbridge", "Vegreville"], data = 0
        )
        self.global_output_df = df(
            index = ["Global_T_Test"],
            columns = ["T_Statistic", "P_Value", "Methylation_Ratio"], data = 0
        )


    def __global_t_test(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame
        ) -> None:
        """
        """
        # Dropping non-data columns
        vegreville_input_df.drop(
            ["#Scaffold", "Bin_Label"], axis = 1, inplace = True
        )
        lethbridge_input_df.drop(
            ["#Scaffold", "Bin_Label"], axis = 1, inplace = True
        )

        # Setting paired T-test results.
        self.global_means_df["Vegreville"] = \
            vegreville_input_df.mean(axis = 0).values
        self.global_means_df["Lethbridge"] = \
            lethbridge_input_df.mean(axis = 0).values
        model = sps.ttest_rel(
            self.global_means_df["Vegreville"],
            self.global_means_df["Lethbridge"]
        )
        self.global_output_df.iloc[0, 0:2] = model[:]

        # Setting methylation ratio.
        self.global_output_df.iloc[0, 2] = vegreville_input_df.sum().sum() / \
            vegreville_input_df.sum().sum()
        self.global_output_df.index = ["global"]


    def global_t_test_and_write(
            self, lethbridge_input_df: pd.DataFrame,
            vegreville_input_df: pd.DataFrame, output_dir_path: str
        ) -> None:
        """
        Perform global paired T-test and save output to a file.
        """
        print("Global t-test start.")

        # Prints stdout to separate file.
        sys.stdout = open("global_t_test_stdout.txt", 'w')
        wrapping_flair = helpers.string_builder((
            '\n', '-' * 5, '*' * 10, '-' * 5, '\n'
        ))
        print(helpers.string_builder((
            wrapping_flair, "Global T-Tests", wrapping_flair
        )))

        self.__set_output_dfs(lethbridge_input_df.columns)
        self.__global_t_test(lethbridge_input_df, vegreville_input_df)
        helpers.write_output(
            output_df = self.global_output_df,
            output_file_name = "global_methylation_ttest.tsv",
            output_dir_path = output_dir_path
        )
        sys.stdout.close()


# Main method.
def paired_t_tests(
        lethbridge_file_path: str, vegreville_file_path: str,
        output_dir_path: str
    ) -> None:
    """
    Performs cross-cultivar and within-cultivar paired T-tests for Lethbridge
    and Vegreville data in parallel using the `multiprocessing` module.
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
    global_output = GlobalPairedTTestOutput()

    # Initialize processes for local and cultivar paired t-test objects.
    func_args = (inputs.lethbridge_df, inputs.vegreville_df, output_dir_path)
    local_process = multiprocessing.Process(
        target = local_output.local_t_test_and_write, args = func_args
    )
    cultivar_process = multiprocessing.Process(
        target = cultivar_output.cultivar_t_test_and_write, args = func_args
    )
    global_process = multiprocessing.Process(
        target = global_output.global_t_test_and_write, args = func_args
    )

    local_process.start()
    cultivar_process.start()
    global_process.start()
    local_process.join()
    cultivar_process.join()
    global_process.join()

    helpers.print_program_runtime("Paired t-test calculations", start_time)


# paired_t_tests(lethbridge_file_path, vegreville_file_path, output_dir_path)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import timeit, df, pd, sps
from . import helpers

import multiprocessing
import sys
import warnings

import statsmodels.formula.api as smf


start_time = timeit.default_timer()
# delta_phenotype_file_path = sys.argv[1]
# delta_methylation_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]


class PhenotypeRegressionInput:
    def __init__(self):
        self.phenotype_df = None
        self.methylation_df = None


    def set_input_dfs(
            self, delta_phenotype_file_path: str,
            delta_methylation_file_path: str
        ) -> None:
        """
        """
        self.phenotype_df = pd.read_table(
            delta_phenotype_file_path, index_col = 0
        )

        cols = ["#Scaffold", "Bin_Label"] + self.phenotype_df.index.tolist()
        self.methylation_df = pd.read_table(
            delta_methylation_file_path, usecols = cols
        )


class PhenotypeRegressionOutput:
    def __init__(self):
        self.current_bin_df = None
        self.phenotype_output_df = None


    def __set_output_dfs(
            self, methylation_input_df: pd.DataFrame
        ) -> None:
        """
        """
        self.current_bin_df = df(
            data = 0, index = methylation_input_df.columns[2:],
            columns = ["delta_phenotype", "delta_methylation"]
        )

        self.phenotype_output_df = methylation_input_df.iloc[:, 0:2]
        self.phenotype_output_df["R_Squared"] = 0
        self.phenotype_output_df["P_Value"] = 0
        self.phenotype_output_df["Significant?"] = False


    def __bin_regression(
            self, phenotype_data: pd.Series, methylation_df: pd.DataFrame
        ) -> None:
        """
        """
        phenotype_label = phenotype_data.name
        self.current_bin_df["delta_phenotype"] = phenotype_data

        for bin_idx in range(methylation_df.shape[0]):
            bin_row = methylation_df.loc[bin_idx]
            quick_search = helpers.string_builder((
                phenotype_label, '-', bin_row[0], '-', str(bin_row[1])
            ))

            methylation_data = pd.to_numeric(bin_row[2:])
            if methylation_data.nonzero()[0].size != 0:
                self.current_bin_df["delta_methylation"] = methylation_data
                warnings.filterwarnings("ignore")
                model = smf.ols(
                    "delta_phenotype ~ delta_methylation",
                    data = self.current_bin_df
                ).fit()

                self.phenotype_output_df.iloc[bin_idx, 2] = model.rsquared
                self.phenotype_output_df.iloc[bin_idx, 3] = 0
                self.phenotype_output_df.iloc[bin_idx, 4] = False
                try:
                    self.phenotype_output_df.iloc[bin_idx, 3] = model.pvalues[1]
                except:
                    pass

                try:
                    self.phenotype_output_df.iloc[bin_idx, 4] = \
                        helpers.significance(model.pvalues.tolist())
                except:
                    pass

                print(
                    helpers.string_builder((
                        '\n', "++++++++++", '\n', quick_search, '\n'
                    ))
                )
                print(model.summary())
                print(helpers.string_builder(('\n', "++++++++++", '\n')))

            else:
                print(
                    helpers.string_builder((
                        '\n', bin_row[0], '-', str(bin_row[1]),
                        " has been filtered for ", phenotype_label, "...", '\n'
                    ))
                )


    def phenotype_regression(
            self, phenotype_data: pd.Series, methylation_input_df: pd.DataFrame,
            output_dir_path: str
        ) -> None:
        """
        """
        phenotype = phenotype_data.name
        sys.stdout = open(
            helpers.string_builder((phenotype, "_stdout.txt")), 'w'
        )
        print(
            helpers.string_builder((
                '\n', "-----**********-----", '\n', "Phenotype: ",
                phenotype, '\n', "-----**********-----", '\n'
            ))
        )

        self.__set_output_dfs(methylation_input_df)
        self.__bin_regression(
            phenotype_data, methylation_input_df
        )

        helpers.write_output(
            self.phenotype_output_df,
            helpers.string_builder((
                phenotype, '_', "phenotype_regression.tsv"
            )), output_dir_path
        )
        sys.stdout.close()


def phenotype_methylation_regression(
        delta_phenotype_file_path: str, delta_methylation_file_path: str,
        output_dir_path: str
    ) -> None:
    """
    Main method.
    """
    print(helpers.string_builder(('\n', "Start.")))
    helpers.remove_trailing_slash([
        delta_phenotype_file_path, delta_methylation_file_path,
        output_dir_path
    ])

    print(helpers.string_builder(('\n', "Setting dataframes...")))
    inputs = PhenotypeRegressionInput()
    inputs.set_input_dfs(delta_phenotype_file_path, delta_methylation_file_path)

    print(helpers.string_builder(('\n', "Performing phenotype regression...")))
    # output = PhenotypeRegressionOutput()
    # output.phenotype_regression(
    #     inputs.phenotype_df["height"], inputs.methylation_df, output_dir_path,
    #     start_time
    # )
    outputs = []
    for phenotype in inputs.phenotype_df.columns.tolist():
        output = PhenotypeRegressionOutput()
        process = multiprocessing.Process(
            target = output.phenotype_regression,
            args = (
                inputs.phenotype_df[phenotype], inputs.methylation_df,
                output_dir_path
            )
        )
        outputs.append(process)

    for process in outputs:
        process.start()

    for process in outputs:
        process.join()

    helpers.print_program_runtime(
        "Phenotype regression analyses", start_time
    )


# phenotype_methylation_regression(
#     delta_methylation_file_path, delta_methylation_file_path, output_dir_path
# )

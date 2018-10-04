#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import timeit, df, pd, sps
from . import helpers

# import sys
import threading

import statsmodels.formula.api as smf


start_time = timeit.default_timer()
# delta_phenotype_file_path = sys.argv[1]
# delta_methylation_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]


class PhenotypeRegressionInput:
    def __init__(self):
        self.phenotype_df = None
        self.methylation_df = None
        self.current_bin_df = None
        self.current_output_df = None


    def __set_dfs(
            self, delta_phenotype_file_path: str,
            delta_methylation_file_path: str
        ) -> None:
        """
        """
        self.phenotype_df = pd.read_table(
            delta_phenotype_file_path, index_col = 0
        )

        cols = ["#Scaffold", "Bin_Label"].append(
            self.phenotype_df.index.tolist()
        )
        self.methylation_df = pd.read_table(
            delta_methylation_file_path, usecols = cols
        )

        self.current_bin_df = df(
            data = 0, index = self.phenotype_df.index,
            columns = ["delta_phenotype", "delta_methylation"]
        )

        self.current_output_df = self.methylation_df.iloc[:, 0:2]
        self.current_output_df["R_Squared"] = 0
        self.current_output_df["P_Value"] = 0
        self.current_output_df["Significant?"] = False


    def __bin_regression(self, phenotype_data: pd.Series) -> None:
        """
        """
        phenotype_label = phenotype_data.name
        self.current_bin_df["delta_phenotype"] = phenotype_data

        for bin_idx in range(self.methylation_df.shape[0]):
            bin_row = self.methylation_df.loc[bin_idx]
            quick_search = helpers.string_builder((
                phenotype_label, "-", bin_row[0], "-", str(bin_row[1])
            ))
            self.current_bin_df["delta_methylation"] = \
                pd.to_numeric(bin_row[2:])

            model = smf.ols(
                "delta_phenotype ~ delta_methylation",
                data = self.current_bin_df
            ).fit()
            self.current_output_df.iloc[bin_idx, 2] = model.rsquared
            self.current_output_df.iloc[bin_idx, 3] = 0
            self.current_output_df.iloc[bin_idx, 4] = False
            try:
                self.current_output_df.iloc[bin_idx, 3] = model.pvalues[1]
            except:
                pass

            try:
                self.current_output_df.iloc[bin_idx, 4] = \
                    helpers.significance(model.pvalues.tolist())
            except:
                pass

            print()
            print("++++++++++")
            print(quick_search)
            print(model.summary())
            print("++++++++++")
            print()


    def __phenotype_regression(self, output_dir_path:str) -> None:
        """
        """
        for phenotype in self.phenotype_df.columns.tolist():
            print()
            print("-----**********-----")
            print(helpers.string_builder(("Phenotype: ", phenotype)))
            print("-----**********-----")
            print()

            self.__bin_regression(self.phenotype_df[phenotype])

            helpers.write_output(
                self.current_output_df,
                helpers.string_builder((
                    phenotype, "_", "phenotype_regression.tsv"
                )), output_dir_path
            )


    def phenotype_methylation_regression(
            self, delta_phenotype_file_path: str,
            delta_methylation_file_path: str, output_dir_path: str
        ) -> None:
        """
        """
        print()
        print("Start")
        helpers.remove_trailing_slash((
            delta_phenotype_file_path, delta_methylation_file_path,
            output_dir_path
        ))

        print()
        print("Setting input dataframes...")
        self.__set_dfs(delta_phenotype_file_path, delta_methylation_file_path)

        print()
        print("Performing phenotype regression...")
        self.__phenotype_regression(output_dir_path)

        helpers.print_program_runtime(
            "Phenotype regression analyses", start_time
        )


# pmr = PhenotypeRegressor()
# pmr.phenotype_methylation_regression(
#     delta_phenotype_file_path, delta_methylation_file_path, output_dir_path
# )

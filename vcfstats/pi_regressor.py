#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple linear regression analysis for DNA and DNAm pi.

Pi_DNAm ~ Pi_DNA + e

The residual becomes the DNAm sequence divergence independent of the effects of
DNA sequence divergence.
"""

from .__init__ import timeit, df, pd
from . import helpers

# import sys

import statsmodels.formula.api as smf


# TODO: move this to user interface.
start_time = timeit.default_timer()
output_df_header = ["Scaff_Bin", "Pi_DNAm", "Pi_DNA", "True_Pi_DNAm"]
# scaffold_list_file_path = sys.argv[1]
# dna_dir_path = sys.argv[2]
# dnam_dir_path = sys.argv[3]
# output_dir_path = sys.argv[4]


class PiRegressor:
    def __init__(self):
        self.bin_dict = {}
        self.scaffold_list = None
        self.current_scaffold_pi_dna = None
        self.current_scaffold_pi_dnam = None
        self.output_df = None


    def __set_current_scaffold(
            self, current_scaffold: str, pi_dna_file_path: str,
            pi_dnam_file_path: str
        ) -> None:
        """
        """
        self.current_scaffold_pi_dna = pd.read_table(pi_dna_file_path)
        self.current_scaffold_pi_dnam = pd.read_table(pi_dnam_file_path)

        for bin_idx in range(self.current_scaffold_pi_dna.shape[0]):
            print(
                helpers.string_builder((
                    "Bin: ", str(self.current_scaffold_pi_dna.loc[bin_idx][0])
                ))
            )
            scaff_bin_lab = helpers.string_builder((
                current_scaffold, "_",
                str(self.current_scaffold_pi_dna.loc[bin_idx][0])
            ))
            pi_dnam = self.current_scaffold_pi_dnam.loc[bin_idx][1]
            pi_dna = self.current_scaffold_pi_dna.loc[bin_idx][1]
            self.bin_dict.__setitem__(scaff_bin_lab, [pi_dnam, pi_dna])


    def __read_scaffold_list(
            self, scaffold_list_file_path: str, dna_dir_path: str,
            dnam_dir_path: str
        ) -> None:
        """
        """
        self.scaffold_list = pd.read_table(
            scaffold_list_file_path, usecols = ["#SCAFFOLD"]
        )

        for idx in range(self.scaffold_list.shape[0]):
            current_scaffold = str(self.scaffold_list.loc[idx][0])
            pi_dna_file_path = helpers.string_builder((
                dna_dir_path, '/', current_scaffold, "_pi_dna.tsv"
            ))
            pi_dnam_file_path = helpers.string_builder((
                dnam_dir_path , '/', current_scaffold, "_pi_dnam.tsv"
            ))

            print(
                helpers.string_builder((
                    "Setting current scaffold bin data row:", '\n',
                    "Scaffold: ", current_scaffold, '\n'
                ))
            )
            self.__set_current_scaffold(
                current_scaffold, pi_dna_file_path, pi_dnam_file_path
            )


    def __set_output_df(self) -> None:
        """
        """
        self.output_df = df(data = self.bin_dict).T.reset_index()
        self.output_df["residuals"] = 0.0
        self.output_df.columns = output_df_header

        # Regression
        model = smf.ols('Pi_DNAm ~ Pi_DNA', data = self.output_df).fit()
        print(helpers.string_builder(('\n', "Model Summary:", '\n')))
        print(model.summary())
        print(model.pvalues)
        self.output_df["True_Pi_DNAm"] = model.resid

        scaff_bin = self.output_df.columns
        self.output_df[["Scaffold", "Bin_Label"]] = \
            self.output_df[scaff_bin[0]].str.split('_', expand = True)

        self.output_df = self.output_df.drop(scaff_bin[0], axis = 1)
        cols = self.output_df.columns.tolist()
        cols = cols[3:] + cols[:3]
        self.output_df = self.output_df[cols]


    def pi_regression(
            self, scaffold_list_file_path: str, dna_dir_path: str,
            dnam_dir_path: str, output_dir_path: str
        ) -> None:
        """
        Main method.
        """
        print(helpers.string_builder(('\n', "Start.")))
        helpers.remove_trailing_slash([
            scaffold_list_file_path, dna_dir_path, dnam_dir_path,
            output_dir_path
        ])

        print(helpers.string_builder(('\n', "Reading scaffold list...")))
        self.__read_scaffold_list(
            scaffold_list_file_path, dna_dir_path, dnam_dir_path
        )

        print(helpers.string_builder(('\n', "Setting output dataframe...")))
        self.__set_output_df()

        helpers.write_output(
            self.output_df, "pi_regression.tsv", output_dir_path
        )

        # Runtime.
        # TODO: move this to user interface as well.
        helpers.print_program_runtime("Pi Regression calculations", start_time)


# pi_rg = PiRegressor()
# pi_rg.pi_regression(
#     scaffold_list_file_path, dna_dir_path, dnam_dir_path, output_dir_path
# )

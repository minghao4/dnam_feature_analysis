#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple regression analysis for DNA and DNAm pi.

Pi_DNAm ~ Pi_DNA + e

The residual becomes the DNAm sequence divergence independent of the effects of
DNA sequence divergence.

"""

import os

from pandas import DataFrame as df
from pandas import read_csv
import statsmodels.formula.api as smf


header_out = ["Scaff_Bin", "Pi_DNAm", "Pi_DNA", "True_Pi_DNAm"]


class PiRegressor:
    def __init__(self):
        self.bin_dict = {}
        self.out_df = df()
        self.curr_pi_dna = df()
        self.curr_pi_dnam = df()
        self.scaff_list = df()


    def __set_curr_scaff(self, curr_scaff, pi_dna_file, pi_dnam_file):
        """

        PiRegressor, String, String, String -> PiRegressor
        """

        self.curr_pi_dna = read_csv(pi_dna_file, sep = '\t')
        self.curr_pi_dnam = read_csv(pi_dnam_file, sep = '\t')

        for bin_idx in range(self.curr_pi_dna.shape[0]):
            print("Bin: " + str(self.curr_pi_dna.loc[bin_idx][0]) + "\n")
            scaff_bin_lab = curr_scaff + "_" + \
                str(self.curr_pi_dna.loc[bin_idx][0])

            pi_dnam = self.curr_pi_dnam.loc[bin_idx][1]
            pi_dna = self.curr_pi_dna.loc[bin_idx][1]

            self.bin_dict.__setitem__(scaff_bin_lab, [pi_dnam, pi_dna])


    def __read_scaff_list(self, scaff_file, input_dir_path):
        """

        PiRegressor, String -> PiRegressor
        """

        self.scaff_list = read_csv(
                              scaff_file,
                              sep = "\t",
                              usecols = ["#SCAFFOLD"]

                          )

        for scaff_idx in range(self.scaff_list.shape[0]):
            curr_scaff = str(self.scaff_list.loc[scaff_idx][0])
            pi_dna_file = input_dir_path + "/" + curr_scaff + "_pi_dna.tsv"
            pi_dnam_file = input_dir_path + "/" + curr_scaff + "_pi_dnam.tsv"

            if not os.path.exists(pi_dna_file):
                continue

            print("Setting current scaffold bin data row:")
            print("Scaffold: " + curr_scaff)
            self.__set_curr_scaff(curr_scaff, pi_dna_file, pi_dnam_file)


    def __set_out_df(self):
        """

        PiRegressor -> PiRegressor
        """

        self.out_df = df(data = self.bin_dict).T.reset_index()
        self.out_df["residuals"] = 0.0
        self.out_df.columns = header_out

        # Regression
        model = smf.ols('Pi_DNAm ~ Pi_DNA', data = self.out_df).fit()
        print(model.summary())
        print(model.pvalues)
        self.out_df["True_Pi_DNAm"] = model.resid

        scaff_bin = self.out_df.columns
        self.out_df[['Scaffold', 'Bin_Label']] = \
            self.out_df[scaff_bin[0]].str.split('_', expand = True)

        self.out_df = self.out_df.drop(scaff_bin[0], axis = 1)
        cols = self.out_df.columns.tolist()
        cols = cols[3:] + cols[:3]
        self.out_df = self.out_df[cols]


    def __write_output(self, output_dir_path):
        """

        PiRegressor, String -> PiRegressor
        """

        output_file = str(output_dir_path + "/" + "pi_regression.tsv")
        if not os.path.isdir(output_dir_path):
            os.mkdir(output_dir_path)

        self.out_df.to_csv(output_file, sep = '\t', index = False)


    def pi_regression(self, scaff_file_path, input_dir_path, output_dir_path):
        """

        PiRegressor, String, String, String -> PiRegressor
        """

        paths = [scaff_file_path, input_dir_path, output_dir_path]
        for path in paths:
            if path.endswith('/'):
                path = path[:-1]

        # if input_dir_path.endswith('/'):
        #     input_dir_path = input_dir_path[:-1]

        # if output_dir_path.endswith('/'):
        #     output_dir_path = output_dir_path[:-1]

        print("Reading scaffold list...\n")
        self.__read_scaff_list(scaff_file_path, input_dir_path)
        self.__set_out_df()
        self.__write_output(output_dir_path)

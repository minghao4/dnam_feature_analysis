#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
"""

from .__init__ import timeit, np, df, pd, sps
from . import helpers

# import sys
from typing import Tuple

start_time = timeit.default_timer()
# bin_file_path = sys.argv[1]
# methylation_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]


class MethylationBinner:
    def __init__(self):
        self.methylation_df = None
        self.bins_output_df = None


    def __set_dfs(
            self, bin_file_path: str, methylation_file_path: str
        ) -> None:
        """
        """
        self.bins_output_df = pd.read_table(bin_file_path)
        self.methylation_df = pd.read_table(methylation_file_path)

        scaffold_position = self.methylation_df.columns[0]
        self.methylation_df[['Scaffold', 'Position']] = \
            self.methylation_df[scaffold_position].str.split('_', expand = True)
        self.methylation_df['Position'] = \
            pd.to_numeric(self.methylation_df['Position'])
        print()
        print("*****")
        print(
            helpers.string_builder((
                "Variants: ", str(self.methylation_df.shape[0])
            ))
        )
        print("*****")
        print()

        # Reordering the columns.
        self.methylation_df = self.methylation_df.drop(
            scaffold_position, axis = 1)
        cols = self.methylation_df.columns.tolist()
        cols = cols[(len(cols) - 2):] + cols[0:(len(cols) - 2)]
        self.methylation_df = self.methylation_df[cols]

        cultivs = cols[2:]
        for cultiv in cultivs:
            self.bins_output_df[cultiv] = 0


    @staticmethod
    def __variant_in_bin(
            current_scaffold: str, bin_lower_bound: int, bin_upper_bound: int,
            variant_scaffold: str, variant_position: int
        ) -> bool:
        """
        """
        in_bin = False
        if variant_scaffold == current_scaffold:
            if variant_position >= bin_lower_bound \
                    and variant_position <= bin_upper_bound:
                in_bin = True

        return in_bin


    def __bin_averaging(
            self, sites: int, current_scaffold: str, bin_idx: int
        ) -> None:
        """
        """
        divide = 1
        if not sites == 0:
            print(
                helpers.string_builder((
                    "Averaging: ", current_scaffold
                ))
            )
            divide = sites

        self.bins_output_df.iloc[bin_idx, 2:] /= divide


    def __read_methylation_df(
            self, bin_idx: int, bookmark: str, current_scaffold: str,
            bin_lower_bound: int, bin_upper_bound: int
        ) -> Tuple[int, int]:
        """
        """
        sites = 0
        for variant_idx in range(bookmark, self.methylation_df.shape[0]):
            variant = self.methylation_df.loc[variant_idx]
            variant_scaffold = variant[0]
            variant_position = variant[1]

            if self.__variant_in_bin(
                    current_scaffold, bin_lower_bound, bin_upper_bound,
                    variant_scaffold, variant_position
                ):
                print(
                    helpers.string_builder((
                        "Adding: ", variant_scaffold, " Position ",
                        str(variant_position)
                    ))
                )
                sites += 1
                bookmark += 1
                self.bins_output_df.iloc[bin_idx, 2:] += variant[2:]

            else:
                self.__bin_averaging(sites, current_scaffold, bin_idx)
                break

        return (bookmark, sites)


    def __read_bin_df(self) -> None:
        """
        """
        methylation_df_bookmark = 0
        sites = 0
        for bin_idx in range(self.bins_output_df.shape[0]):
            current_bin = self.bins_output_df.loc[bin_idx]
            current_scaffold = current_bin[0]
            current_bin_label = float(current_bin[1])
            bin_lower_bound = 0
            bin_upper_bound = 0

            print()
            print(
                helpers.string_builder((
                    "Reading: ", current_scaffold, " Bin ",
                    str(current_bin_label)
                ))
            )

            # TODO: fix this, 200bp is laziness because default bin is 400
            if current_bin_label % 200 == 0:
                bin_lower_bound = current_bin_label - 200 + 1
                bin_upper_bound = current_bin_label + 200

            else:
                bin_lower_bound = \
                    current_bin_label - current_bin_label % 200 + 1
                bin_upper_bound = \
                    current_bin_label + current_bin_label % 200 + 1

            bookmark_sites = self.__read_methylation_df(
                bin_idx, methylation_df_bookmark, current_scaffold,
                bin_lower_bound, bin_upper_bound
            )
            methylation_df_bookmark = bookmark_sites[0]
            sites = bookmark_sites[1]
            print()
            print("-----")
            print(methylation_df_bookmark)
            print("-----")
            print()

            if methylation_df_bookmark == self.methylation_df.shape[0]:
                break

        self.__bin_averaging(sites, current_scaffold, bin_idx)


    def calculate_all_bin_methylation(
            self, bin_file_path: str, methylation_file_path: str,
            output_dir_path: str
        ) -> None:
        """
        Main method.
        """
        print()
        print("Start.")
        helpers.remove_trailing_slash([
            bin_file_path, methylation_file_path, output_dir_path
        ])

        print()
        print("Setting input dataframes...")
        self.__set_dfs(bin_file_path, methylation_file_path)
        # sys.exit()

        print("Calculating average methylation...")
        print()
        self.__read_bin_df()


        helpers.write_output(
            self.bins_output_df, "methylation_bins.tsv", output_dir_path
        )

        helpers.print_program_runtime("Methylation binning", start_time)


# mb = MethylationBinner()
# mb.calculate_all_bin_methylation(
#     bin_file_path, methylation_file_path, output_dir_path
# )

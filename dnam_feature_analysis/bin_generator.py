#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Objective: generate a sorted bins file from given sorted scaffold sizes file.

Inputs:
- Sorted scaffold name and scaffold sizes TSV file path.
- Output directory path.

Outputs:
- TSV file holding sorted bins file.

Bin size 400bp as shortest scaffolds are 400bp long.

"""

from . import List, timeit, df, math, np, pd
from . import helpers

file_header = ["#Scaffold", "Bin_Label"]


class BinGenerator:
    def __init__(self) -> None:
        self.scaffold_df = None
        self.bin_df = None


    def __set_dfs(
            self, header: List[str], scaffold_sizes_file_path: str
        ) -> None:
        """
        Sets the input and output dataframes.
        """
        self.scaffold_df = pd.read_table(scaffold_sizes_file_path)
        self.bin_df = df(data = 0, columns = header)


    def __iter_scaffolds(self, row: pd.Series, header: List[str]) -> None:
        """
        Given a scaffold, create bins and append to the bin dataframe.
        """
        # Bin size 400bp
        scaffold_name = row[0]
        scaffold_size = int(row[1])
        num_bins = math.ceil(scaffold_size / 400)

        # Final bin
        final_bin_label = int((scaffold_size % 400) / 2) # Initialize
        if final_bin_label != 0: # If final bin is incomplete
            final_bin_label += 400 * (num_bins - 1)

         # If final bin is 1 (e.g. scaffold length is 401bp)
        elif scaffold_size % 400 == 1:
            final_bin_label == scaffold_size

        # Bin labels are midpoints of the bin
        scaffold_bins = df(
            (np.arange(num_bins * 2).reshape(num_bins, 2) * 200).astype(float),
            columns = header
        )
        scaffold_bins["#Scaffold"] = scaffold_name
        scaffold_bins[num_bins - 1, 1] = final_bin_label
        self.bin_df.append(scaffold_bins, ignore_index = True)
        del scaffold_bins # mem management


    # Main method.
    def bin_generator(
            self, scaffold_sizes_file_path: str,
            output_dir_path: str
        ) -> None:
        """
        Generate 400bp bins.
        """
        start_time = timeit.default_timer()
        scaffold_sizes_file_path, output_dir_path = \
            helpers.remove_trailing_slash((
                scaffold_sizes_file_path, output_dir_path
            ))

        print("\nStart.\nSetting dataframes...")
        self.__set_dfs(file_header, scaffold_sizes_file_path)

        print("\nGenerating bins...")
        tmp = self.scaffold_df.apply(
            self.__iter_scaffolds, args = (file_header,), axis = 1
        )
        del tmp

        helpers.write_output(self.bin_df, "sorted_bins.tsv", output_dir_path)
        helpers.print_program_runtime("Bin generation", start_time)

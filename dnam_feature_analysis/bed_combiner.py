#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Objective: to combine BED files from all cultivars.

Inputs:
- List of cultivar names
- Paths to Lethbridge and Vegreville cultivar BED file directories

Output:
- Write combined output for each location as a TSV file

Output Columns:
- Scaffold_Position
- 1 column per cultivar, holding the beta value (methylation level) at the
  genomic location.

"""

from . import multiprocessing, List, sys, timeit, Tuple, natsorted, df, pd
from . import helpers


class BedCombiner:
    def __init__(self, location_label: str, bed_dir_path: str) -> None:
        self.location_label = location_label
        self.bed_dir_path = helpers.remove_trailing_slash((bed_dir_path))
        self.cultivar_df = None
        self.output_df = df()


    def __read_current_cultivar_file(
            self, cultivar: str, cultivar_bed_file_path: str
        ) -> None:
        """
        Reads the current cultivar BED file.
        """
        self.cultivar_df = pd.read_table(
            cultivar_bed_file_path, names = ["#Scaffold", "Position", cultivar],
            usecols = [0, 2, 7] # Scaffold, position, and beta value
        )


    def __index_cultivar_df(self) -> None:
        """
        Combines the "#Scaffold" and "Position" columns and sets the resulting
        column as the index of the current cultivar dataframe.
        """
        self.cultivar_df.index = self.cultivar_df["#Scaffold"] + '_' + \
            self.cultivar_df["Position"].map(str)

        self.cultivar_df.drop(
            ["#Scaffold", "Position"], axis = 1, inplace = True
        )

    # `natsort` module used for reindexing.
    def __concat_cultivar_output_dfs(self) -> None:
        """
        Concatenates the current cultivar dataframe to the output dataframe,
        reindexing the output dataframe.
        """
        self.output_df = pd.concat([self.output_df, self.cultivar_df], axis = 1)

        print("Reindexing output dataframe...")
        self.output_df = self.output_df.reindex(
            index = natsorted(self.output_df.index)
        )


    def loc_bed_combiner(self, cultivars: List[str]) -> None:
        """
        Performs the steps for combining all cultivar BED files at the current
        location iteratively.
        """
        print(helpers.string_builder((
            self.location_label, " BED combining start."
        )))

        # Prints stdout to a separate file.
        sys.stdout = open(
            helpers.string_builder((
                self.location_label, "_bed_combine_stdout.txt"
            )), 'w'
        )

        print("Looping through cultivar BED files...")
        for cultivar in cultivars:
            cultivar_bed_file_path = helpers.string_builder((
                self.bed_dir_path, "/cultivars/", cultivar, '_',
                self.location_label, ".bed"
            ))

            print(helpers.string_builder(("\nCurrently reading: ", cultivar)))
            self.__read_current_cultivar_file(cultivar, cultivar_bed_file_path)

            print(helpers.string_builder(("Arranging index for ", cultivar)))
            self.__index_cultivar_df()

            print(helpers.string_builder((
                "Concatenating ", cultivar, " data to output dataframe..."
            )))
            self.__concat_cultivar_output_dfs()

            helpers.write_output(
                output_df = self.output_df,
                output_file_name = "sorted_methylation_levels.tsv",
                output_dir_path = self.bed_dir_path, write_index = True
            )

        sys.stdout.close()


# Main method.
def bed_combiner(self, cultivars: List[str], bed_dir_paths: Tuple[str]) -> None:
    """
    Combines BED files at Lethbridge and Vegreville in parallel
    using the `multiprocessing` module.
    """
    start_time = timeit.default_timer() # Initialize starting time.

    print("\nStart\n") # Initialize BedCombiner objects for the two locations.
    lethbridge_bed_combiner = BedCombiner('L', bed_dir_paths[0])
    vegreville_bed_combiner = BedCombiner('V', bed_dir_paths[1])

    # Initialize a process for each location's BedCombiner object.
    lethbridge_process = multiprocessing.Process(
        target = lethbridge_bed_combiner.loc_bed_combiner, args = (cultivars,)
    )
    vegreville_process = multiprocessing.Process(
        target = vegreville_bed_combiner.loc_bed_combiner, args = (cultivars,)
    )

    # Start processes and rejoin them when complete.
    lethbridge_process.start()
    vegreville_process.start()
    lethbridge_process.join()
    vegreville_process.join()

    helpers.print_program_runtime("Combining BED files ", start_time)

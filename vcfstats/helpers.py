#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Helper functions.
"""

from .__init__ import os, timeit, pd

from typing import List, Tuple


def significance(model: list) -> bool:
    """
    """
    significance = False
    if model[0] != 0 and model[1] != 0:
        if model[1] <= 0.05:
            significance = True

    return significance


def string_builder(str_list: Tuple[str]) -> str:
    return "".join(str_list)


def remove_trailing_slash(file_paths_list: Tuple[str]) -> Tuple[str]:
    """
    """
    for file_path in file_paths_list:
        if file_path.endswith("/"):
            file_path = file_path[:-1]

    return file_paths_list


def create_output_directory(output_dir_path: str) -> None:
    """
    """
    if not os.path.isdir(output_dir_path):
        print("Creating output directory...")
        os.mkdir(output_dir_path)


def write_output(
        output_df: pd.DataFrame, output_file_name: str, output_dir_path: str,
        write_index: bool = False
    ) -> None:
    """
    """
    output_file = string_builder((output_dir_path, '/', output_file_name))
    create_output_directory(output_dir_path)
    print()
    print(string_builder(("Writing ", output_file, " to ", output_dir_path)))
    output_df.to_csv(output_file, sep = '\t', index = write_index)


def print_program_runtime(program_name: str, start_time: float) -> None:
    """
    """
    raw_runtime = timeit.default_timer() - start_time
    runtime = int(raw_runtime)
    hours = 0
    minutes = 0
    seconds = runtime % 60
    if runtime >= 3600:
        hours = int(runtime / 3600)
        minutes = int(runtime % 3600 / 60)

    elif runtime >= 60:
        minutes = int(runtime / 60)

    # Looks like this:
    #
    # ==========
    # (Output text)
    # ==========
    #
    wrapping_flair = string_builder(('\n', '=' * 10, '\n'))
    print(string_builder((
            wrapping_flair, program_name, " complete.\n", "Raw runtime: ",
            str(raw_runtime), "s.\n", "Program runtime: ", str(hours), "h ",
            str(minutes), "m ", str(seconds), "s.\n", wrapping_flair
    )))

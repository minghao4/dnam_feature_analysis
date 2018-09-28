#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Helper functions.
"""

import os
import timeit
from typing import Tuple, List

import pandas as pd


def string_builder(str_list: Tuple[str]) -> str:
    return "".join(str_list)


def remove_trailing_slash(file_paths_list: List[str]) -> List[str]:
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
        df: pd.DataFrame, output_file_name: str, output_dir_path: str
    ) -> None:
    """
    """
    output_file = string_builder((output_dir_path, '/', output_file_name))
    create_output_directory(output_dir_path)
    print()
    print(string_builder(("Writing ", output_file, " to ", output_dir_path)))
    df.to_csv(output_file, sep = '\t', index = False)


def print_program_runtime(program_name: str, start_time: float) -> None:
    """
    """
    raw_runtime = timeit.default_timer() - start_time
    runtime = int(raw_runtime)
    hours = 0
    minutes = 0
    seconds = runtime % 60
    if runtime >= 3600:
        hours = runtime / 3600
        minutes = runtime % 3600 / 60

    elif runtime >= 60:
        minutes = runtime / 60

    print()
    print("==========")
    print(string_builder((program_name, " complete.")))
    print(string_builder(("Raw runtime: ", str(raw_runtime), "s.")))
    print(
        string_builder((
            "Program runtime: ", str(hours), "h ", str(minutes), "m ",
            str(seconds), "s."
        ))
    )
    print("==========")
    print()

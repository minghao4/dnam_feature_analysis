#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import timeit

"""
Helpers
"""

def remove_trailing_slash(file_paths_list):
    """

    List[String] -> List[String]
    """

    for file_path in file_paths_list:
        if file_path.endswith("/"):
            file_path = file_path[:-1]

    return file_paths_list


def create_output_directory(output_dir_path):
    """

    String -> None
    """

    if not os.path.isdir(output_dir_path):
        os.mkdir(output_dir_path)


def print_program_runtime(program_name, start_time):
    """

    String, Float -> None
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
    print(program_name + " complete.")
    print("Raw runtime: " + str(raw_runtime) + "s.")
    print("Program runtime: " + \
        str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s.")
    print("==========")
    print()


def write_output(output_file_name, output_dir_path):
    """

    String, String -> None
    """

    output_file = str(output_dir_path + "/" + output_file_name)
    self.out_
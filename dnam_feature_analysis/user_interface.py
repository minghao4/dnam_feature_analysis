#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
User interface.

"""

import argparse
from . import bed_combiner, bin_generator, delta_methylation_and_phenotype, \
    helpers, methylation_binner, paired_t_tester, phenotype_regressor

cultivars = [
    "canda", "cfx1", "cfx2", "crs1", "delores", "finola", "grandi",
    "joey", "katani", "picolo", "silesia", "x59"
]


def user_interface():
    program_description = \
        "DNA Methylation Feature Analysis of Lethbridge and Vegreville Plants"
    parser = argparse.ArgumentParser(description = program_description)

    # Add arguments
    bed_combiner_help = "Combine BED files."
    parser.add_argument(
        "-bc", "--bed_combiner", type = str, nargs = 2,
        metavar = ("lethbridge_directory", "vegreville_directory"),
        default = None, help = bed_combiner_help
    )

    bin_generator_help = "Generate bins."
    parser.add_argument(
        "-bg", "--bin_generator", type = str, nargs = 2,
        metavar = ("scaffold_sizes_file", "output_directory"),
        default = None, help = bin_generator_help
    )

    methylation_binner_help = "Calculate bin methylation."
    parser.add_argument(
        "-mb", "--methylation_binner", type = str, nargs = 3,
        metavar = ("sorted_bins_file", "methylation_file", "output_directory"),
        default = None, help = methylation_binner_help
    )

    paired_t_tester_help = "Paired t-tests."
    parser.add_argument(
        "-ptt", "--paired_t_tester", nargs = 3,
        metavar = ("lethbridge_file", "vegreville_file", "output_directory"),
        default = None, help = paired_t_tester_help
    )

    delta_mp_help = "Delta - Vegreville minus Lethbridge"
    parser.add_argument(
        "-dmp", "--delta_mp", type = str, nargs = 5,
        metavar = (
            "lethbridge_methylation_file", "vegreville_methylation_file",
            "lethbridge_phenotype_file", "vegreville_phenotype_file",
            "output_directory"
        ), default = None, help = delta_mp_help
    )

    phenotype_regressor_help = "Phenotype regression."
    parser.add_argument(
        "-pr", "--phenotype_regressor", type = str, nargs = 3,
        metavar = (
            "delta_phenotype_file", "delta_methylation_file", "output_directory"
        ), default = None, help = phenotype_regressor_help
    )

    args = parser.parse_args()

    # Arguments
    if args.bed_combiner != None:
        bed_combiner.bed_combiner(cultivars, args[0], args[1])

    elif args.bin_generator != None:
        bg_obj = bin_generator.BinGenerator()
        bg_obj.bin_generator(args[0], args[1])

    elif args.methylation_binner != None:
        mb_obj = methylation_binner.MethylationBinner()
        mb_obj.calculate_all_bin_methylation(args[0], args[1], args[2])

    elif args.paired_t_tester != None:
        paired_t_tester.paired_t_tests(args[0], args[1], args[2])

    elif args.delta_mp != None:
        delta_methylation_and_phenotype.delta(
            args[0], args[1], args[2], args[3], args[4]
        )

    elif args.phenotype_regressor != None:
        phenotype_regressor.phenotype_methylation_regression(
            args[0], args[1], args[2]
        )

    else:
        parser.print_help()

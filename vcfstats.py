#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
# import vcfstats.vcf_parser as vp
# import vcfstats.sequence_divergence as sd
# import vcfstats.pi_regression as pr
# import vcfstats.methylation_bin as mb
from vcfstats import global_t_test as gtt

# input_path = sys.argv[1]
# parser = vp.VcfParser()

# parser.parse_all_vcfs(input_path)

# bin_width = int(sys.argv[1])
# var_file = sys.argv[2]
# scaff_sizes_file = sys.argv[3]
# output_dir_path = sys.argv[4]
# header_out = ["#Distance", "Pi"]

# pi_calc = sd.PiCalculator()
# pi_calc.calculate_pi_all_scaffolds(
#     bin_width,
#     var_file,
#     scaff_sizes_file,
#     header_out,
#     output_dir_path,

# )

# scaff_file_path = sys.argv[1]
# dna_dir_path = sys.argv[2]
# dnam_dir_path = sys.argv[3]
# output_dir_path = sys.argv[4]

# pi_rg = pr.PiRegressor()
# pi_rg.pi_regression(
#     scaff_file_path,
#     dna_dir_path,
#     dnam_dir_path,
#     output_dir_path,

# )

# bin_file = sys.argv[1]
# lethbridge_methLvl_file = sys.argv[2]
# vegreville_methLvl_file = sys.argv[3]
# lethbridge_output_folder = sys.argv[4]
# vegreville_output_folder = sys.argv[5]

# m_bin = mb.MethylationBinner()
# m_bin.calculate_all_bin_methylation(bin_file, lethbridge_methLvl_file, lethbridge_output_folder)
# m_bin.calculate_all_bin_methylation(bin_file, vegreville_methLvl_file, vegreville_output_folder)


lethbridge_methLvl_file = sys.argv[1]
vegreville_methLvl_file = sys.argv[2]
output_folder = sys.argv[3]

g_tt = gtt.GlobalTTester()
g_tt.global_t_testing(lethbridge_methLvl_file, vegreville_methLvl_file, output_folder)


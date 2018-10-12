#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
# import vcfstats.vcf_parser as vp
# import vcfstats.sequence_divergence as sd
# import vcfstats.pi_regression as pr
# import vcfstats.methylation_binner as mb
# from vcfstats import global_t_test as gtt
from vcfstats import cultivar_t_tester as ctt
# from vcfstats import phenotype_regressor as pmr

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

# bin_file_path = sys.argv[1]
# methylation_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]

# m_bin = mb.MethylationBinner()
# m_bin.calculate_all_bin_methylation(bin_file_path, methylation_file_path, output_dir_path)


# lethbridge_methLvl_file = sys.argv[1]
# vegreville_methLvl_file = sys.argv[2]
# output_folder = sys.argv[3]

# g_tt = gtt.GlobalTTester()
# g_tt.global_t_testing(lethbridge_methLvl_file, vegreville_methLvl_file, output_folder)

lethbridge_file_path = sys.argv[1]
vegreville_file_path = sys.argv[2]
output_dir_path = sys.argv[3]

ctt.paired_t_tests(lethbridge_file_path, vegreville_file_path, output_dir_path)

# delta_phenotype_file_path = sys.argv[1]
# delta_methylation_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]

# pmr.phenotype_methylation_regression(
#     delta_phenotype_file_path, delta_methylation_file_path, output_dir_path
# )

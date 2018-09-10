#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import vcfstats.vcf_parser as vp
# import vcfstats.sequence_divergence as sd
# import vcfstats.pi_regression as pr

input_path = sys.argv[1]
parser = vp.VcfParser()

parser.parse_all_vcfs(input_path)

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
# input_dir_path = sys.argv[2]
# output_dir_path = sys.argv[3]

# pi_rg = pr.PiRegressor()
# pi_rg.pi_regression(scaff_file_path, input_dir_path, output_dir_path)

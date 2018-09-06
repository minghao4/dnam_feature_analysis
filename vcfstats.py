#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import vcfstats.vcf_parser as vp

input_path = sys.argv[1]
parser = vp.VcfParser()

parser.parse_all_vcfs(input_path)

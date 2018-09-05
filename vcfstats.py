import sys
from vcfstats.vcf_parser import VcfParser

input_path = sys.argv[1]
parser = VcfParser()

parser.parse_all_vcfs(input_path)

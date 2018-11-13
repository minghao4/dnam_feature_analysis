# DNA Methylation Feature Analysis

## Usage

``` python

python dnam_feature_analysis/
usage:  [-h] [-bc lethbridge_directory vegreville_directory]
        [-bg scaffold_sizes_file output_directory]
        [-mb sorted_bins_file methylation_file output_directory]
        [-ptt lethbridge_file vegreville_file output_directory]
        [-dmp lethbridge_methylation_file vegreville_methylation_file lethbridge_phenotype_file vegreville_phenotype_file output_directory]
        [-pr delta_phenotype_file delta_methylation_file output_directory]

DNA Methylation Feature Analysis of Lethbridge and Vegreville Plants

optional arguments:
  -h, --help            show this help message and exit
  -bc lethbridge_directory vegreville_directory, --bed_combiner lethbridge_directory vegreville_directory
                        Combine BED files.
  -bg scaffold_sizes_file output_directory, --bin_generator scaffold_sizes_file output_directory
                        Generate bins.
  -mb sorted_bins_file methylation_file output_directory, --methylation_binner sorted_bins_file methylation_file output_directory
                        Calculate bin methylation.
  -ptt lethbridge_file vegreville_file output_directory, --paired_t_tester lethbridge_file vegreville_file output_directory
                        Paired t-tests.
  -dmp lethbridge_methylation_file vegreville_methylation_file lethbridge_phenotype_file vegreville_phenotype_file output_directory, --delta_mp lethbridge_methylation_file vegreville_methylation_file lethbridge_phenotype_file vegreville_phenotype_file output_directory
                        Delta - Vegreville minus Lethbridge
  -pr delta_phenotype_file delta_methylation_file output_directory, --phenotype_regressor delta_phenotype_file delta_methylation_file output_directory
                        Phenotype regression.

```

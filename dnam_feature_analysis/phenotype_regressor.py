#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Objective: perform simple linear regression where delta phenotype is regressed
against delta methylation at each bin.

Inputs:
- Delta phenotype TSV file path.
- Delta methylation TSV file path.
- Output directory path.

Outputs:
- Log TXT stdout files.
- TSV file holding regression results (R Squared value, p-value, nominal
  significance) for each phenotype.

"""

from . import multiprocessing, sys, timeit, warnings, df, pd, smf, sps
from . import helpers


# delta_phenotype_file_path = sys.argv[1]
# delta_methylation_file_path = sys.argv[2]
# output_dir_path = sys.argv[3]


class PhenotypeRegressionInput:
    def __init__(self) -> None:
        self.phenotype_df = None
        self.methylation_df = None


    def set_input_dfs(
            self, delta_phenotype_file_path: str,
            delta_methylation_file_path: str
        ) -> None:
        """
        Set input dataframes.
        """
        self.phenotype_df = pd.read_table(
            delta_phenotype_file_path, index_col = 0
        )

        cols = ["#Scaffold", "Bin_Label"] + self.phenotype_df.index.tolist()
        self.methylation_df = pd.read_table(
            delta_methylation_file_path, usecols = cols
        )


class PhenotypeRegressionOutput:
    def __init__(self) -> None:
        self.phenotype_output_df = None


    def __set_output_df(self, methylation_input_df: pd.DataFrame) -> None:
        """
        Set output dataframe.
        """
        # Default values.
        self.phenotype_output_df = methylation_input_df.iloc[:, 0:2]
        self.phenotype_output_df["R_Squared"] = 0
        self.phenotype_output_df["P_Value"] = 0
        self.phenotype_output_df["Significant?"] = False


    def __iter_bin_regression(
            self, bin_row: pd.Series, phenotype_data: pd.Series,
            phenotype_label: str, bin_df_index: pd.Index
        ) -> None:
        """
        Perform simple linear regression on delta methylation and delta
        phenotype.
        """
        bin_idx = bin_row.name
        quick_search = helpers.string_builder((
            phenotype_label, '-', bin_row[0], '-', str(bin_row[1])
        ))
        methylation_data = pd.to_numeric(bin_row[2:])

        if methylation_data.nonzero()[0].size != 0:
            # Initialize bin dataframe.
            current_bin_df = df(
                data = 0, index = bin_df_index,
                columns = ["delta_phenotype", "delta_methylation"]
            )
            current_bin_df["delta_phenotype"] = phenotype_data
            current_bin_df["delta_methylation"] = methylation_data

            # Perform simple linear regression (ordinary least squares)
            warnings.filterwarnings("ignore")
            model = smf.ols(
                "delta_phenotype ~ delta_methylation",
                data = current_bin_df
            ).fit()
            del current_bin_df # mem management

            self.phenotype_output_df.iloc[bin_idx, 2] = model.rsquared
            self.phenotype_output_df.iloc[bin_idx, 3] = 0
            self.phenotype_output_df.iloc[bin_idx, 4] = False
            try:
                self.phenotype_output_df.iloc[bin_idx, 3] = model.pvalues[1]
            except:
                pass

            try:
                self.phenotype_output_df.iloc[bin_idx, 4] = \
                    helpers.significance(model.pvalues.tolist())
            except:
                pass

            wrapping_flair = helpers.string_builder(('\n', '+' * 10, '\n'))
            print(helpers.string_builder((
                wrapping_flair, quick_search, '\n'
            )))
            print(model.summary())
            print(wrapping_flair)

        else:
            print(helpers.string_builder((
                '\n', bin_row[0], '-', str(bin_row[1]),
                " has been filtered for ", phenotype_label, "...", '\n'
            )))


    def __bin_regression(
            self, phenotype_data: pd.Series, methylation_df: pd.DataFrame
        ) -> None:
        """
        Perform simple linear regression on delta methylation and delta
        phenotype for all bins.
        """
        phenotype_label = phenotype_data.name
        bin_df_index = methylation_df.columns[2:]
        tmp = methylation_df.apply(
            self.__iter_bin_regression,
            args = (phenotype_data, phenotype_label, bin_df_index), axis = 1
        )
        del tmp


    def phenotype_regression(
            self, phenotype_data: pd.Series, methylation_input_df: pd.DataFrame,
            output_dir_path: str
        ) -> None:
        """
        Perform simple linear regression on delta methylation and delta
        phenotype for the current phenotype.
        """
        phenotype = phenotype_data.name
        print(helpers.string_builder((phenotype, "start.")))

        # Prints stdout to separate file.
        sys.stdout = open(
            helpers.string_builder((phenotype, "_stdout.txt")), 'w'
        )
        wrapping_flair = helpers.string_builder((
            '\n', '-' * 5, '*' * 10, '-' * 5, '\n'
        ))
        print(helpers.string_builder((
            wrapping_flair, "Phenotype: ", phenotype, wrapping_flair
        )))

        self.__set_output_df(methylation_input_df)
        self.__bin_regression(phenotype_data, methylation_input_df)
        helpers.write_output(
            self.phenotype_output_df,
            helpers.string_builder((
                phenotype, '_', "phenotype_regression.tsv"
            )), output_dir_path
        )
        sys.stdout.close()


# Main method.
def phenotype_methylation_regression(
        delta_phenotype_file_path: str, delta_methylation_file_path: str,
        output_dir_path: str
    ) -> None:
    """
    Perform simple linear regression on delta methylation and delta
    phenotype for all phenotypes within the delta phenotype file.
    """
    start_time = timeit.default_timer()
    delta_phenotype_file_path, delta_methylation_file_path, output_dir_path = \
        helpers.remove_trailing_slash((
            delta_phenotype_file_path, delta_methylation_file_path,
            output_dir_path
        ))

    print("\nStart.\nSetting dataframes...") # Initialize input object.
    inputs = PhenotypeRegressionInput()
    inputs.set_input_dfs(delta_phenotype_file_path, delta_methylation_file_path)

    print("\nPerforming phenotype regression...") # Initialize output objects.
    output_obs = {}
    output_processes = {} # Initialize output processes.
    for phenotype in inputs.phenotype_df.columns.tolist():
        output_obs[phenotype] = PhenotypeRegressionOutput()
        process = multiprocessing.Process(
            target = output_obs[phenotype].phenotype_regression,
            args = (
                inputs.phenotype_df[phenotype], inputs.methylation_df,
                output_dir_path
            )
        )
        output_processes[phenotype] = process

    for phenotype in output_processes:
        output_processes[phenotype].start()

    for phenotype in output_processes:
        output_processes[phenotype].join()

    helpers.print_program_runtime("Phenotype regression analyses", start_time)


# phenotype_methylation_regression(
#     delta_phenotype_file_path, delta_methylation_file_path, output_dir_path
# )

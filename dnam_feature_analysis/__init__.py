#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DNA Methylation Feature Analysis.

"""

__all__ = [
    "bed_combiner", "bin_generator", "delta_methylation_and_phenotype",
    "helpers", "methylation_binner", "paired_t_tester", "phenotype_regressor",
    "user_interface"
]

# Native python libs
import math
import multiprocessing
import os
import sys
import timeit
from typing import List, Tuple
import warnings

# External libs
from natsort import natsorted
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import scipy.stats as sps
import statsmodels.formula.api as smf

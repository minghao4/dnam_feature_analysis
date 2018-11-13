#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VCF Stats

"""

import multiprocessing
import os
import sys
import timeit
from typing import List, Tuple
import warnings
import statsmodels.formula.api as smf

import math
from natsort import natsorted
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import scipy.stats as sps
import statsmodels.formula.api as smf

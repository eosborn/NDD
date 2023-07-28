#goal: summarize frequency of truncating mutations in NVA cohort

#input: genes_counts.pkl
#output: dataframe, rows = genes,
    #columns = (gene, total_count, num_samples, fs_del, fs_insertion, startloss, stopgain, stoploss)

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import pickle
import os
import sys

from truncating_summary import *

if __name__ == "__main__":
    #parse command line input
    #format:
    #python truncating_summary.py genes_counts.pkl
    input = pd.read_pickle(sys.argv[1])

    output = input.loc[:, ['gene', 'name', 'cross_reference',
                           'fs_insertion_count','fs_del_count',
                           'startloss_count', 'stopgain_count', 'stoploss_count']]

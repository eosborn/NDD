#plotting count data

from plotnine import ggplot, geom_histogram, aes, scale_x_continuous, scale_y_continuous, labs, ggtitle
import pandas as pd
import numpy as np
import pyarrow.feather as feather
import sys
import os

#parse command line input
    #format:
    #python plotting.py <in_directory> <out_directory>
#in directory contains feather files of count data of grouped variant files
in_directory = sys.argv[1]
out_directory = sys.argv[2]

for file_name in os.listdir(in_directory):
    if file_name.endswith("SUB.feather"):
        file_path = in_directory + file_name
        df = feather.read_feather(file_path)

        freqs = pd.DataFrame(df["count"].value_counts())
        freqs = freqs.sort_index()
        
        #names for figure
        name = file_name.split(".")[0]
        item = file_name.split("_")[0]
        if item[-1] == "s":
            item = item[:-1]
        title = item + " frequency in *more than one sample*"
        x_title = "number of unique " + item + "s that have a given frequency" 
        y_title = "frequency of " + item + " across samples"

        my_plot = (ggplot(df, aes(x = "count"))
                + geom_histogram(bins = freqs.index[-1], color="darkblue", fill = "lightblue")
                + scale_x_continuous(breaks = freqs.index)
                + scale_y_continuous(breaks = freqs["count"])
                + ggtitle(title)
                + labs(y = y_title, x = x_title))
        out_file =  out_directory + name + "_frequency_hist.png"
        my_plot.save(out_file)


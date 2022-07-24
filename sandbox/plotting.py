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
in_directory = sys.argv[1]
out_directory = sys.argv[2]

for file_name in os.listdir(in_directory):
    if file_name.endswith(".feather"):
        name = file_name.split("_")[0]
        file_path = in_directory + file_name
        df = feather.read_feather(file_path)

        freqs = pd.DataFrame(df["count"].value_counts())
        freqs = freqs.sort_index()

        my_plot = (ggplot(df, aes(x = "count"))
                + geom_histogram(bins = freqs.index[-1], color="darkblue", fill = "lightblue")
                + scale_x_continuous(breaks = freqs.index)
                + scale_y_continuous(breaks = freqs["count"])
                + ggtitle("frequency of " + name + " counts")
                + labs(y = "frequency", x = "number of times a given "+ name +" is counted"))
        out_file =  out_directory + name + "_frequency_hist.png"
        my_plot.save(out_file)

##        #print the genes that have the higher counts
##        gene_list = []
##        for count in freqs.index:
##            if count != 1:
##                gene_index = np.where(df["count"]==count)[0][0]
##                gene_list.append(df.iloc[gene_index]["Gene_name"])
##        print(gene_list)

##df_subset = df[df["count"] != 1]
##df_out = df_subset.to_csv(sep='\t', index=False)
##writefile = open("test_data/counts/Genes_multiple_hits.txt", "w")
##writefile.write(df_out)
##writefile.close

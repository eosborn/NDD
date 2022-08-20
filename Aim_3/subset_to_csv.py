#goal: create csv ouptut of gene/grouping data
    #removing most common genes (TTN, MUC12, FLG in my case)

import pandas as pd
import numpy as np
import sys
import os
import pickle

from subset_to_csv import *

if __name__ == "__main__":
    #parse command line input
    #format:
    #python subset_to_csv.py <in_directory> <out_directory>
    #in_directory will contain all the count files as .pickle files
    #out_directory will be where csv output files will be written
    #genes_rm = list of genes to remove
    in_directory = sys.argv[1]
    out_directory = sys.argv[2]
    #genes_rm = input("Input list of genes to remove: ")
    genes_rm = ['TTN', 'MUC12', 'FLG']

    #create a dictionary of count files
    print("\n", "Reading in count files", "\n")
    d = create_dict(in_directory)

    #subsetting count dfs to remove most common genes
    print("removing {} from count dataframes".format(genes_rm), "\n")
    #d_2 is a dictionary where output dfs will be stored 
    d_2 = {}
    for df_name in d:
        #special case for gene df --> just remove gene rows, write out to csv 
        if df_name == "genes":
            #probably a better way to this 
            d_2[df_name] =  d[df_name][~d[df_name]['genes'].isin(genes_rm)].reset_index()
            d_2[df_name] = d_2[df_name].drop(columns = ["index"])
        #grouping count dfs --> need to loop through dictionaries in 'samples' columns,
            # remove genes from value lists ("N") 
            # subract the number of genes removed from the overall 'count' value ('count' - N)
        else:
            d_2[df_name] = remove_genes(d[df_name], genes_rm)

    #writing out csvs
    print("writing csv output files", "\n")
    for df in d_2:
        return_df = d_2[df].to_csv(index=False)
        out_file_name = out_directory + df + "_subset.csv"
        try:
            writefile = open(out_file_name, "w")
        except FileNotFoundError:
            os.mkdir(out_directory)
            writefile = open(out_file_name, "w")
        writefile.write(return_df)
        writefile.close
            

########## functions ##########
    
def create_dict(in_directory):
    """
    """
    d = {}
    for file in os.listdir(in_directory):
        if file.endswith(".pkl"):
            print(file)
            file_name = file.split("_")[0]
            d[file_name] = pd.read_pickle(in_directory + file)
    return(d)

def remove_genes(count_df, genes_rm):
    """
    """
    for index, row in count_df.iterrows():
        #count = row[1]
        #samples_dict = row[2] {"sample_name" : ["gene1", "gene2", etc]}
        # N will count the number of times the removed genes occur per row
        N = 0 
        for sample in row[2]:
            #count occurences of genes to remove and subract that from count
            for gene in genes_rm:
                N += row[2][sample].count(gene)
        row[1] -= N
        #then remove genes from samples_dict 
        for sample in row[2]:
            try:
                row[2][sample] = [elem for elem in row[2][sample] if elem not in genes_rm]
            except ValueError:
                pass
    return(count_df)

#goal: create files that count genes and groupings across all filtered and grouped variant files

#v2 changes:
#1. pickle output instead of feather because better for python-to-python data transfer
#2. using __name__ variable to make script more readable/linear
#3. created additional functions to simplify/streamline main code
#4. only output groupings with more than one observation across samples
#5. added try/except clause to create output directory if needed
#6. added metadata columns to output
#7. sort output by descending count 

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import pickle
import os
import sys

from counting_v2 import *

if __name__ == "__main__":
    #parse command line input
    #format:
    #python counting_v2.py <in_directory> <out_directory> <meta_data_directory>
    #in_directory will contain all the grouped variant files as .feather files
    #meta_data_directory will contain gene/grouping meta data as tsv files,
        #file names formated as "only_*_tsv.txt"
    in_directory = sys.argv[1]
    out_directory = sys.argv[2]
    meta_directory = sys.argv[3]

    #create a dictionary of dfs that counts will be added to
    things_to_count = ["genes", "pathways", "networks", "modules", "brite"]
    d = {}
    for item in things_to_count:
        d[item] = create_empty_df(item)

    #count each sample in input directory
    #example file_name: "18PA008847_ASD-genes.2022-07-05.pDet.grouped.feather"
    print("Looping through variant files and counting", "\n")
    for file_name in os.listdir(in_directory):
        #only want to read in feather files 
        if file_name.endswith(".feather"):
            sample_df = read_in_df(file_name, in_directory)
            sample_name = file_name.split("_")[0]
            #for each sample file, count the five categories and store data in "d" df dictionary 
            for count_file in d:
                #special case to count the genes 
                if count_file == "genes":
                    #don't need to explode because each row represents a single gene/variant observation
                    d[count_file] = count_genes(sample_df["Gene_name"], d[count_file], sample_name)
                #to count groupings 
                else:
                    #explode grouping so that each row represents a single grouping observation
                    df_to_count = sample_df[["Gene_name", count_file]].explode(count_file)
                    d[count_file] = count_groupings(df_to_count, d[count_file], count_file, sample_name)
        else:
            pass

    #create dictionary of metadata
    #first ask for where metadata tsv files are located
    print("Loading in metadata", "\n")
    #meta_dir = str(input("Input metadata tsv directory: "))
    meta_d = {}
    for file_name in os.listdir(meta_directory):
        #only want to read in tsv files 
        if file_name.endswith("tsv.txt"):
            #file name example: "only_brite_tsv.txt"
            #print(file_name)
            #item_name = str(input("Input item name (ie 'genes', 'pathways', etc: "))
            item_name = file_name.split("_")[1]
            file_path = meta_directory + file_name
            meta_d[item_name] = create_meta_df(item_name, file_path)
        else:
            pass

    #remove groupings that appear in just one sample
    print("Performing final editing of output files", "\n")
    for count_file in d:
        tmp_file = remove_single_counts(d[count_file])
        tmp_file = remove_single_samples(tmp_file, count_file)
        #add in metadata
        for item in meta_d:
            if item == count_file:
                d[count_file] = add_metadata(tmp_file, meta_d[item])
            else:
                pass
        #sort by counts
        out_file = tmp_file.sort_values(by=["count"], ascending=False, ignore_index = True)

        #write out to txt file and pickle
        return_df = out_file.to_csv(sep='\t', index=False)
        out_file_name = out_directory + count_file + "_counts.txt"
        try:
            writefile = open(out_file_name, "w")
        except FileNotFoundError:
            os.mkdir(out_directory)
            writefile = open(out_file_name, "w")
        writefile.write(return_df)
        writefile.close
        out_file.to_pickle(out_directory + count_file + "_counts.pkl")

                                  
########## functions ##########
    
def create_empty_df(name):
    """
    """
    empty_df = pd.DataFrame(columns = [name, "count", "samples"])
    return(empty_df)

def read_in_df(file_name, in_directory):
    """
    """
    file_path = in_directory + file_name
    df = feather.read_feather(file_path)
    #rename gene column to remove the period
    df = df.rename(columns={"Gene.refGene" : "Gene_name"})
    return(df)

def count_genes(df_to_count, count_df, sample_name):
    """
    """
    #get list of items that already have rows in the count df 
    counted_already = list(count_df["genes"])

    #loop through genes
    for item in df_to_count:
        if item not in counted_already:
            tmp_df = pd.DataFrame([[item , 1 , {sample_name :1 }]],
                                  columns = ["genes", "count", "samples"])
            #used .concat instead of .append because latter is deprecated
            count_df = pd.concat([count_df, tmp_df], ignore_index=True)
            #add newly counted item to counted_already list
            counted_already.append(item)
        else:
            count_df.loc[count_df["genes"] == item, "count"] += 1
            index = int(np.where(count_df["genes"]==item)[0][0])
            #if this is the first time this sample has been counted for this particular gene 
            if sample_name not in count_df.iloc[index]["samples"]:
                #add sample to samples dictionary
                count_df.iloc[index]["samples"][sample_name] = 1
            #if this is not the first time 
            else:
                #add 1 to samples dictionary
                count_df.iloc[index]["samples"][sample_name] += 1
    return(count_df)

def count_groupings(df_to_count, count_df, item_being_counted, sample_name):
    """
    """
    #get list of items that already have rows in the count df 
    counted_already = list(count_df[item_being_counted])

    #iterate over df_to_count rows so that the gene_names and groupings are kept together
    for index, row in df_to_count.iterrows():
        #skip over rows with no groupings 
        if type(row[item_being_counted]) != float:
            if row[item_being_counted] not in counted_already:
                #add new row to count_df
                tmp_df = pd.DataFrame([[row[item_being_counted] , 1 , {sample_name : [row["Gene_name"]]}]],
                                  columns = [item_being_counted, "count", "samples"])
                count_df = pd.concat([count_df, tmp_df], ignore_index=True)
                #add newly counted item to counted_already list
                counted_already.append(row[item_being_counted])
            else:
                #add 1 to grouping count
                count_df.loc[count_df[item_being_counted] == row[item_being_counted], "count"] += 1
                index_count_df = int(np.where(count_df[item_being_counted]==row[item_being_counted])[0][0])

                #if this is the first time this grouping has appeared in this sample --> add new key:value
                if sample_name not in count_df.iloc[index_count_df]["samples"]:
                    #add sample to samples dictionary with gene 
                    count_df.iloc[index_count_df]["samples"][sample_name] = [row["Gene_name"]]
                #if this is not the first time --> extend sample value list  
                else:
                    #add gene name to sample value list 
                    count_df.iloc[index_count_df]["samples"][sample_name].append(row["Gene_name"])
        else:
            pass
    return(count_df)

def remove_single_counts(count_df):
    """
    """
    file = count_df.loc[count_df["count"] != 1]
    file = file.reset_index()
    return(file)

def remove_single_samples(tmp_df, name_of_item):
    """
    """
    #create new df that will hold appended lines from tmp_df
    out_df = create_empty_df(name_of_item)
    for index, row in tmp_df.iterrows():
        if len(tmp_df.iloc[index]["samples"]) > 1:
            row_to_append = pd.DataFrame([[row[name_of_item] , row["count"] , row["samples"]]],
                                      columns = [name_of_item, "count", "samples"])
            out_df = pd.concat([out_df, row_to_append], ignore_index=True)
    return(out_df)

def create_meta_df(item_name, file_path):
    """
    """
    if item_name.lower() == "genes":
        df = pd.read_csv(file_path, sep = '\t', \
                         names = ["symbol", "KEGG_entry_ID", "name"], dtype=str)
    else:
        df = pd.read_csv(file_path, sep = '\t', \
                         names = ["KEGG_entry_ID", "name", "num_NDD_genes"], dtype=str)
    return(df)
   
def add_metadata(output_df, meta_df):
    """
    """
    #new cols for genes: 'KEGG_entry_ID', 'name'
    #new cols for groupings: 'name', 'num_NDD_genes'
    
    output_df[meta_df.columns[1]] = ""
    output_df[meta_df.columns[2]] = ""
    for index, row in output_df.iterrows():
        try:
            meta_df_index = np.where(meta_df.iloc[0:,0]==row[0])[0][0]
            output_df.iloc[index,3] = meta_df.iloc[meta_df_index, 1]
            output_df.iloc[index,4] = meta_df.iloc[meta_df_index, 2]
        except IndexError:
            output_df.iloc[index,3] = str(input("Input {} {} KEGG entry ID: "\
                                                .format(output_df.columns[0],row[0])))
            output_df.iloc[index,4] = str(input("Input {} {} name: "\
                                                .format(output_df.columns[0],row[0])))
    return(output_df)

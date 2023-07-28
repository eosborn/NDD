#goal: create files that count genes and groupings across all filtered and grouped variant files

#v3 changes:
#1. reordering of first four columns: gene/<group>, KEGG_entry_ID, name, total_count
#2. new column for gene count df: Xref.refGene (unique to gene)
#3. additional 11 new columns for gene count df:
    # one for each ExonicFunc.refGene annotation (frameshift insertion, stopgain, etc)
#4. each sample has its own column (<sample_ID>_count)

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import pickle
import os
import sys

from counting_v3 import *

if __name__ == "__main__":
    #parse command line input
    #format:
    #python counting_v2.py <in_directory> <meta_data_directory> <out_directory>
    #in_directory will contain all the grouped variant files as .feather files
    #meta_data_directory on NDD located at ~/Aim_2a/parsed_to_tsv
        #contain gene/grouping meta data as tsv files
        #file names formated as "only_*_tsv.txt"
    in_directory = sys.argv[1]
    meta_directory = sys.argv[2]
    out_directory = sys.argv[3]

    #create a dictionary of dfs that counts will be added to
    things_to_count = ["genes", "pathways", "networks", "modules", "brite"]
    d = {}
    for item in things_to_count:
        d[item] = create_empty_df(item)

    #create dictionary of metadata
    print("Loading in metadata", "\n")
    meta_d = {}
    for file_name in os.listdir(meta_directory):
        #only want to read in tsv files 
        if file_name.endswith("tsv.txt"):
            #file name example: "only_brite_tsv.txt"
            #item_name = str(input("Input item name (ie 'genes', 'pathways', etc: "))
            item_name = file_name.split("_")[1]
            file_path = meta_directory + file_name
            meta_d[item_name] = create_meta_df(item_name, file_path)
        else:
            pass

    #count each sample in input directory
    #example file_name: "18PA008847_ASD-genes.2022-07-05.pDet.grouped.feather"
    print("Looping through variant files and counting", "\n")
    for file_name in os.listdir(in_directory):
        #only want to read in feather files 
        if file_name.endswith(".feather"):
            sample_df = read_in_df(file_name, in_directory)
            sample_name = file_name.split("_")[0]
            print("\n", sample_name)
            #for each sample file, count the five categories and store data in "d" df dictionary 
            for count_file in d:
                #special case to count the genes 
                if count_file == "genes":
                    #don't need to explode because each row represents a single gene/variant observation
                    d[count_file] = count_genes(sample_df, d[count_file], sample_name, meta_d[count_file])
                #to count groupings 
                else:
                    d[count_file] = count_groupings(sample_df, d[count_file], count_file, sample_name, meta_d[count_file])
                    #(old approach)explode grouping so that each row represents a single grouping observation
                    #df_to_count = sample_df[["Gene_name", count_file]].explode(count_file)
                    #d[count_file] = count_groupings(df_to_count, d[count_file], count_file, sample_name, meta_d[count_file])
        else:
            pass

    #remove groupings that appear in just one sample
    print("Performing final editing of output files", "\n")
    for count_file in d:
        tmp_file = remove_single_counts(d[count_file])
        #tmp_file = remove_single_samples(tmp_file, count_file)
        #sort by counts
        out_file = tmp_file.sort_values(by=["total_count"], ascending=False, ignore_index = True).reset_index()

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
    if name == "genes":
        empty_df = pd.DataFrame(columns = ["gene", "KEGG_entry_ID", "name", "total_count",
                                           "cross_reference",
                                           "fs_insertion_count", "fs_del_count",
                                           "startloss_count", "stopgain_count", "stoploss_count",
                                           "nonfs_insertion_count", "nonfs_del_count", "nonfs_block_sub_count",
                                           "nonsyn_SNV_count", "syn_SNV_count", "unknown_count"])
    else:
        empty_df = pd.DataFrame(columns = ["KEGG_entry_ID", "name", "total_count", "genes"])
    return(empty_df)

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

def read_in_df(file_name, in_directory):
    """
    """
    file_path = in_directory + file_name
    df = feather.read_feather(file_path)
    #rename input columns to remove the period
    df = df.rename(columns={"Gene.refGene" : "Gene_name",
                            "ExonicFunc.refGene" : "variant_function",
                            "Xref.refGene" : "cross_reference"})
    return(df)

def count_genes(df_to_count, count_df, sample_name, meta_df):
    """
    """
    #add sample name as column to count_df
    count_df[sample_name] = 0
    
    #get list of items that already have rows in the count df 
    counted_already = list(count_df["gene"])

    #loop through sample file that needs to be counted 
    for index,row in df_to_count.iterrows():
        if row.Gene_name not in counted_already:
            #create empty df row
            tmp_df = pd.DataFrame(np.nan, index = [0], columns = count_df.columns.tolist())
            #fill in columns in empty row
            tmp_df.gene = row.Gene_name
            try:
                tmp_df.KEGG_entry_ID = meta_df[meta_df["symbol"] == row.Gene_name].KEGG_entry_ID.item()
                tmp_df.name = meta_df[meta_df["symbol"] == row.Gene_name].name.item()
            except ValueError:
                tmp_df.KEGG_entry_ID = ""
                tmp_df.name = ""
                #tmp_df.KEGG_entry_ID = str(input("Input {} KEGG entry ID: "\
                #                                .format(row.Gene_name)))
                #tmp_df.name = str(input("Input {} name: "\
                #                                .format(row.Gene_name)))
                #add this info to the meta_df so don't have to repeat
                tmp_meta_df = pd.DataFrame([[row.Gene_name, tmp_df.KEGG_entry_ID.item(), tmp_df.name.item()]],
                                           columns = ["symbol", "KEGG_entry_ID", "name"])
                meta_df = pd.concat([meta_df, tmp_meta_df], ignore_index=True)
            tmp_df.total_count = 1
            tmp_df.cross_reference = row.cross_reference
            #0s for variant_function columns
            tmp_df.iloc[:,5:16] = int(0)
            #0s for all other samples except the current one
            tmp_df.iloc[:,16:-1] = int(0)
            tmp_df[sample_name] = 1
            #update count for the type of variant
            if row.variant_function == "frameshift insertion":
                tmp_df.fs_insertion_count = 1
            elif row.variant_function == "frameshift deletion":
                tmp_df.fs_del_count = 1
            elif row.variant_function == "startloss":
                tmp_df.startloss_count = 1
            elif row.variant_function == "stopgain":
                tmp_df.stopgain_count = 1
            elif row.variant_function == "stoploss":
                tmp_df.stoploss_count = 1
            elif row.variant_function == "nonframeshift deletion":
                tmp_df.nonfs_del_count = 1
            elif row.variant_function == "nonframeshift insertion":
                tmp_df.nonfs_insertion_count = 1
            elif row.variant_function == "nonframeshift block substitution":
                tmp_df.nonfs_block_sub_count = 1
            elif row.variant_function == "nonsynonymous SNV":
                tmp_df.nonsyn_SNV_count = 1
            elif row.variant_function == "synonymous SNV":
                tmp_df.syn_SNV_count = 1
            elif row.variant_function == "unknown":
                tmp_df.unknown_count = 1
            else:
                print("unrecorded variant:", sample_name, row.Gene_name, row.variant_function)
            #used .concat instead of .append because latter is deprecated
            count_df = pd.concat([count_df, tmp_df], ignore_index=True)
            #add newly counted item to counted_already list
            counted_already.append(row.Gene_name)
        else:
            #add 1 to total_count column and to sample column
            count_df.loc[count_df["gene"] == row.Gene_name, "total_count"] += 1
            count_df.loc[count_df["gene"] == row.Gene_name, sample_name] += 1
            #update variant function count
            if row.variant_function == "frameshift insertion":
                count_df.loc[count_df["gene"] == row.Gene_name, "fs_insertion_count"] += 1
            elif row.variant_function == "frameshift deletion":
                count_df.loc[count_df["gene"] == row.Gene_name, "fs_del_count"] += 1
            elif row.variant_function == "startloss":
                count_df.loc[count_df["gene"] == row.Gene_name, "startloss_count"] += 1
            elif row.variant_function == "stopgain":
                count_df.loc[count_df["gene"] == row.Gene_name, "stopgain_count"] += 1
            elif row.variant_function == "stoploss":
                count_df.loc[count_df["gene"] == row.Gene_name, "stoploss_count"] += 1
            elif row.variant_function == "nonframeshift deletion":
                count_df.loc[count_df["gene"] == row.Gene_name, "nonfs_del_count"] += 1
            elif row.variant_function == "nonframeshift insertion":
                count_df.loc[count_df["gene"] == row.Gene_name, "nonfs_insertion_count"] += 1
            elif row.variant_function == "nonframeshift block substitution":
                count_df.loc[count_df["gene"] == row.Gene_name, "nonfs_block_sub_count"] += 1
            elif row.variant_function == "nonsynonymous SNV":
                count_df.loc[count_df["gene"] == row.Gene_name, "nonsyn_SNV_count"] += 1
            elif row.variant_function == "synonymous SNV":
                count_df.loc[count_df["gene"] == row.Gene_name, "syn_SNV_count"] += 1
            elif row.variant_function == "unknown":
                count_df.loc[count_df["gene"] == row.Gene_name, "unknown_count"] += 1
            else:
                print("unrecorded variant:", sample_name, row.Gene_name, row.variant_function)
            # (old code, might need later) index = int(np.where(count_df["gene"]==item)[0][0])
    return(count_df)

def count_groupings(df_to_count, count_df, item_being_counted, sample_name, meta_df):
    """
    """
    #add sample name as column to count_df
    count_df[sample_name] = 0

    #get list of items that already have rows in the count df 
    counted_already = list(count_df["KEGG_entry_ID"])

    #iterate over df_to_count rows so that the gene_names and groupings are kept together
    for index, row in df_to_count.iterrows():
        grouping_list = row[item_being_counted]
        if len(grouping_list) > 0:
            #check each group ID contained in the list
            for group_ID in grouping_list:            
                if group_ID not in counted_already:
                    #create empty row
                    tmp_df = pd.DataFrame(np.nan, index = [0], columns = count_df.columns.tolist())
                    tmp_df.KEGG_entry_ID = group_ID
                    tmp_df.name = meta_df[meta_df["KEGG_entry_ID"] == group_ID].name.item()
                    tmp_df.total_count = 1
                    tmp_df["genes"] = pd.Series([[row.Gene_name]], dtype=object)
                    #0s for all previous sample counts, add one for current sample
                    tmp_df.iloc[:,4:-1] = int(0)
                    tmp_df[sample_name] = 1
                    #add new row to count_df
                    count_df = pd.concat([count_df, tmp_df], ignore_index=True)
                    #add newly counted item to counted_already list
                    counted_already.append(group_ID)
                else:
                    #add 1 to grouping count and sample count
                    count_df.loc[count_df["KEGG_entry_ID"] == group_ID, "total_count"] += 1
                    count_df.loc[count_df["KEGG_entry_ID"] == group_ID, sample_name] += 1
                    #check to see if gene name is already in genes list, add it if it is not
                    if row.Gene_name in count_df.loc[count_df["KEGG_entry_ID"] == group_ID, "genes"].item():
                        pass
                    else:
                        condition = count_df["KEGG_entry_ID"] == group_ID
                        index = count_df[condition].index.item()
                        #print("\n", count_df.at[index, "genes"])
                        count_df.at[index, "genes"].append(row.Gene_name)
        else:
            pass
    return(count_df)

def remove_single_counts(count_df):
    """
    """
    file = count_df.loc[count_df["total_count"] != 1]
    file = file.reset_index()
    return(file)

##def remove_single_samples(tmp_df, name_of_item):
##    """
##    """
##    #create new df that will hold appended lines from tmp_df
##    out_df = create_empty_df(name_of_item)
##    for index, row in tmp_df.iterrows():
##        if len(tmp_df.iloc[index]["samples"]) > 1:
##            row_to_append = pd.DataFrame([[row[name_of_item] , row["count"] , row["samples"]]],
##                                      columns = [name_of_item, "count", "samples"])
##            out_df = pd.concat([out_df, row_to_append], ignore_index=True)
##    return(out_df)

#goal: create files that count genes and groupings across all filtered and grouped variant files

import pandas as pd
import numpy as np
import pyarrow.feather as feather
import os
import sys


def create_empty_dict(name):
    """
    """
    empty_df = pd.DataFrame(columns = [name, "count", "samples"])
    return(empty_df)

def count_genes(df_to_count, count_df, sample_name):
    """
    """
    #get list of items that already have rows in the count df 
    counted_already = list(count_df["Gene_name"])

    #loop through genes
    for item in df_to_count:
        if item not in counted_already:
            tmp_df = pd.DataFrame([[item , 1 , {sample_name :1 }]],
                                  columns = ["Gene_name", "count", "samples"])
            #used .concat instead of .append because latter is deprecated
            count_df = pd.concat([count_df, tmp_df], ignore_index=True)
            #add newly counted item to counted_already list
            counted_already.append(item)
        else:
            count_df.loc[count_df["Gene_name"] == item, "count"] += 1
            index = int(np.where(count_df["Gene_name"]==item)[0][0])
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


##########
#parse command line input
    #format:
    #python counting.py <in_directory> <out_directory>
#in directory will contain all the grouped variant files as .feather files 
in_directory = sys.argv[1]
out_directory = sys.argv[2]

#create a dictionary of dfs that counts will be added to 
things_to_count = ["Gene_name", "pathways", "networks", "modules", "brite"]
d = {}
for item in things_to_count:
    d[item] = create_empty_dict(item)

for file_name in os.listdir(in_directory):
    #only want to read in feather files 
    if file_name.endswith(".feather"):
        sample_name = file_name.split(".")[0]
        file_path = in_directory + file_name
        df = feather.read_feather(file_path)
        #rename gene column to remove the period
        df = df.rename(columns={"Gene.refGene" : "Gene_name"})
        
        #for each sample file, count the five categories and store data in d dictionary 
        for count_file in d:
            #special case to count the genes 
            if count_file == "Gene_name":
                #don't need to explode because each row represents a single gene/variant observation
                df_to_count = df.Gene_name
                d[count_file] = count_genes(df_to_count, d[count_file], sample_name)
            #to count groupings 
            else:
                #explode grouping so that each row represents a single grouping observation
                df_to_count = df[["Gene_name", count_file]].explode(count_file)
                d[count_file] = count_groupings(df_to_count, d[count_file], count_file, sample_name)
    else:
        pass


#remove groupings that appear in just one sample
d_2 = {}
for count_file in d:
    file = d[count_file]
    file = file.loc[file["count"] != 1]
    file = file.reset_index()

    d_2[count_file] = create_empty_dict(count_file)
    for index, row in file.iterrows():
        if len(file.iloc[index]["samples"]) > 1:
            tmp_df = pd.DataFrame([[row[count_file] , row["count"] , row["samples"]]],
                                  columns = [count_file, "count", "samples"])
            d_2[count_file] = pd.concat([d_2[count_file], tmp_df], ignore_index=True)
    

#os.mkdir(out_directory)
for count_file in d:
    return_df = d[count_file].to_csv(sep='\t', index=False)
    out_file_name = out_directory + count_file + "_counts.txt"
    writefile = open(out_file_name, "w")
    writefile.write(return_df)
    writefile.close

    feather.write_feather(d[count_file], out_directory + count_file + "_counts.feather")

for count_file in d_2:
    return_df = d_2[count_file].to_csv(sep='\t', index=False)
    out_file_name = out_directory + count_file + "_countsSUB.txt"
    writefile = open(out_file_name, "w")
    writefile.write(return_df)
    writefile.close

    feather.write_feather(d_2[count_file], out_directory + count_file + "_countsSUB.feather")
    

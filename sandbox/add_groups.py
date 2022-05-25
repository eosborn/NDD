import pandas as pd
import numpy as np
import json

#load in sample filtered variant file
df = pd.read_csv("18PA004616.ASD-genes.2022-05-05.pDet.txt", sep='\t')
#print(df.head())

#create empy new column that in a future step will store grouping data
#df["pathways"] = ""
#print(df.head())

#function that takes as input a variant and outputs grouping list taken from gene:[groupings] dict
#maybe eventually want to generalize to get_groups
def get_pathways(variant):
    """ Return the list of groups associated with the input varaint's gene

        Parameters
        ----------
            variant: <class 'pandas.core.series.Series'>

        Return
        __________
            variant_list : list
                The list of grouping names associated with the input variant's gene name
            [] : list
                An empty list for a gene that does not exist in the gene:[grouping] dictionary
    """
    dict_file = open("../KEGG_output/dicts/gene_dicts/pathways_gene_dict.txt")
    dict_string = dict_file.readline()
    dict_file.close()

    dict_string = dict_string.replace("\'","\"")
    dict = json.loads(dict_string)
    #error --> Expecting property name enclosed in double quotes: line 1 column 2 (char 1)
    #issue with quotations --> JSON only allows enclosing string with double quotes
    # (before) "{'ADA':...
    # (after) '{"ADA":...

    while True:
        try:
            variant_list = dict[variant[6]]
            return(variant_list)
        except KeyError:
            return([])
            
#if __name__ == '__main__':
#    print(get_pathways.__doc__)


#apply get_pathways function to each row of df
print("Original DataFrame:\n", df.head())
df["pathways"] = df.apply(lambda row : get_pathways(row), axis = 1)
print("New DataFrame:\n", df.head())
    
    

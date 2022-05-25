import pandas as pd
import numpy as np
import json

def add_groups(sample_variant_file, **kwargs):
    """ Return the input sample variant file with new columns for each of the input groupings

        Parameters
        ----------
            sample_variant_file : <class 'str'>
                A tab delimited text file containing variants for a single sample/genome
            **kwargs :
                Keyword argument = name of grouping (pathways, networks, etc)
                Variable = text file containing dictionary {key = gene name ; value = list of associated groupings}
        Return
        ----------
            df : <class 'pandas.core.frame.DataFrame'>
                The inputed sample variant df will contain a new column for each input grouping
                Each new column will contain the lists of associated group names for the particular gene/variant given in each row
    """
    df = pd.read_csv(sample_variant_file, sep='\t')

    d = {}
    for key, value in kwargs.items():
        dict_file = open(value)
        dict_string = dict_file.readline()
        dict_file.close()

        dict_string = dict_string.replace("\'","\"")
        dict = json.loads(dict_string)
        d[key] = dict

    #now have a pandas df of variants, and a dictionary of dictionaries {key = grouping_name : value = grouping dictionary}
    #next, want to loop through each dictionary, adding list of group names to each row of the 
    for group in d:
        dict = d[group]
        add_columns(df,group,dict)

    return_df = df.to_csv(sep='\t')
    out_file_name = sample_variant_file[:-3] + "grouped.txt"
    writefile = open(out_file_name, "w")
    writefile.write(return_df)
    writefile.close
    return(df.head())   

def add_columns(df,group,dict):
    """ return the input df with a new column for the input grouping name and dict

        Parameters
        ----------
            df : <class 'pandas.core.frame.DataFrame'>
            group : <class 'str'>
            d[group] : <class 'dict'>
    
        Return
        ----------
            df : <class 'pandas.core.frame.DataFrame'>
                The input df with an new column named <group> and containing lists of grouping names 
    """
    df[group] = df.apply(lambda row : add_group_lists(row,dict), axis = 1)
    return(df)

def add_group_lists(variant, dict):
    """ Return the list of groups associated with the input varaint's gene

        Parameters
        ----------
            variant : <class 'pandas.core.series.Series'>
            dict : <class 'dict'>

        Return
        ----------
            groups_list : <class 'list'>
                The list of grouping names associated with the input variant's gene name
            [] : <class 'list'>
                An empty list for a gene that does not exist in the gene:[grouping] dictionary
    """
    while True:
        try:
            groups_list = dict[variant[6]]
            return(groups_list)
        except KeyError:
            return([])
            

##if __name__ == '__main__':
##    print(add_groups.__doc__)
##    print(add_columns.__doc__)
##    print(add_group_lists.__doc__)
    

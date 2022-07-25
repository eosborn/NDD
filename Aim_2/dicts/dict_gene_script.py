#import idlelib.idle

#goal: create dictionaries for four parsed KEGG output files
#key = gene : value = [pathways]

#working directory:
import os
os.chdir('/Users/Ellen/Files/Shankar/NDD/pathways/KEGG_output/parsed')
#os.getcwd()

#loop through files in directory
directory = '/Users/Ellen/Files/Shankar/NDD/pathways/KEGG_output/parsed'

for file in os.listdir(directory):
    filename = str(file)
    file = open(file)
    line = file.readline()
    line = line[:-1]
    dict = {}
    while line != "":
        if line[:3] == "hsa" or line[0] == "0":
            pathway_name = ""
            for character in line:
                if character != " ":
                    pathway_name += character
                else:
                    break
            line = file.readline()
            line = line[:-1]
            while line != "" and (line[:3] != "hsa" and line[0] != "0"):
                if line not in dict:
                    dict[line] = [pathway_name]
                else:
                    dict[line].append(pathway_name)
                line = file.readline()
                line = line[:-1]
    file.close()
    dict_file_output = filename[14:-11]+"_gene_dict.txt"
    print(dict_file_output)
    print(len(dict))
    writefile = open("/Users/Ellen/Files/Shankar/NDD/pathways/KEGG_output/dicts/gene_dicts/{}".format(dict_file_output),"w")
    writefile.write(str(dict))
    writefile.close()
    #print(dict)

#check --> grep each gene name (dict keys) in parsed txt file, see if matches len of that gene's value pathway list 

    

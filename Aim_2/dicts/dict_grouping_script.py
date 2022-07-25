#import idlelib.idle

#goal: create dictionaries for four parsed KEGG output files
#current working directory:
import os
os.chdir('/Users/Ellen/Files/Shankar/NDD/pathways/KEGG_output/parsed')
#os.getcwd()

#formatting
#brite --> hsa03036 (43 classifications)
#modules --> hsa_M00065 (91 modules)
#networks --> 06460 (all networks start with 0) (115 networks)
#pathways --> hsa0110 (337 pathways)

#loop through files in directory
directory = '/Users/Ellen/Files/Shankar/NDD/pathways/KEGG_output/parsed'

for file in os.listdir(directory):
    if file[-3:] == "txt":
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
                gene_list = []
                line = file.readline()
                while line != "" and (line[:3] != "hsa" and line[0] != "0"):
                    line = line[:-1]
                    gene_list.append(line)
                    line = file.readline()
                dict[pathway_name] = gene_list
        #print(dict)
        dict_name = str(filename[14:-11])+"_pathway_dict.txt"
##        print(filename)
##        print(len(dict))
        file.close()
        writefile = open("/Users/Ellen/Files/Shankar/NDD/pathways/KEGG_output/dicts/{}".format(dict_name),"w")
        writefile.write(str(dict))
        writefile.close()

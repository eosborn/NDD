### Aim 2: scrape KEGG mapper output for pathways/modules and create dictionaires to organize groupings

#### *dicts/*
.txt files that contain dictionaires of scraped and parsed KEGG output
- ##### *gene_dicts/*
gene names are keys and associated grouping lists are values; created with `dict_gene_script.py` located in parent directory
	- these files were chosen for downstream analysis because it was easier to lookup groupings associated with a given gene while looping through variant files 
- ##### *grouping dicts/*
grouping names are keys and associated gene lists are values; created with `dict_grouping_script.py` 
	- these files were not used in downstream analysis 

#### *parsed/*
.txt files contain grouping ID and name followed by associated genes taken from gene list generated in Aim 1
- parsed from `NDD_Gene_List_*_.txt` files located in parent directory

#### *parsed_to_tsv/*
.txt files contain more gene/grouping information than the *parsed/* .txt files, tab-delimited for downstream python use
- created by `create_tsv.sh` and parsed from `NDD_Gene_List_*_.txt` files located in parent directory
- `genes_tsv` columns = gene symbol, gene ID, and gene name 
- `only_*_tsv` columns = grouping ID, grouping name, and number of genes from gene list associated with the given grouping
- ##### *col_files/* intermediary files pasted together to create tsv files 
- ##### *grouping_genes_files/* intermediary files of all genes associated with each grouping, used to create `genes_tsv` in parent directory
	- there is overlap in genes found in each grouping class which was removed to create `genes_tsv`
- ##### `create_tsv.sh` used to create tsv files 
running `create_tsv.sh`: 
1. enter `bash create_tsv.sh` in command line/terminal
2. enter input files at prompt (used `NDD_Gene_List_*.txt` files), 
3. enter output name that user would like to use for each given input file at prompt (ie 'brite', 'pathway', etc)
*output automatically placed in parsed_to_tsv/ and script will automatically relocate to parsed_to_tsv/ when done running*   

#### `NDD_Gene_List_*.txt` files 
scraped output of [KEGG Mapper search](https://www.genome.jp/kegg/mapper/search.html) 

To do:
1. make python scripts executable from the command line with user input of parsed files
2. make python scripts output pickle files to make downstream implementations easier 
3. add bash script used to scrape KEGG output and create parsed txt files


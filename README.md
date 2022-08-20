# Variant and Gene Identification in Neurodevelopmental Delay Cohort 
## project for the [Translational Health Data Science Fellowship](https://datalab.ucdavis.edu/2021/12/17/announcing-2022-translational-health-data-science-fellows/) through UC Davis' DataLab

### Aim 1: collect list of genes associated with neurodevelopemntal disorders (ASD, speech delay, nonverbal)
- see `NDD_Gene_List.txt` in *Aim_1/* directory
- generated using a broad literature search

To do:
1. add literature sources  
	
### Aim 2: scrape KEGG mapper output for pathways/modules and create dictionaires to organize groupings
- raw KEGG mapper output located in *Aim_2/* as .txt files
- parsed output (just grouping names and gene names) located in *Aim_2/parsed/*
- dictionaries with gene names as *keys* located in *Aim_2/dicts/gene_dicts/*
	- made using `dict_gene_script.py` python script
- dictionaries with grouping names as *keys* located in *Aim_2/dicts/grouping_dicts* 
	- made using `dict_grouping_script.py` python script
- tab deliminted files containing metadata on genes and groups located in *Aim_2/parsed_to_tsv/*
	- made using `create_tsv.sh`   	

To do:
1. include shell script used to parse KEGG output prior to dictionary creation

### Aim 3: apply gene list and groupings to existing sequencing data for patient cohort
- see `add_groups.py` and `counting_v2.py` python scripts located in *Aim_3/* directory

To do:
1. add *README* to directory 
2. fill in function docstrings in `counting.py`

### Aim 4: identify relevant genes/groupings for clinical interpretation  
- currently workshopping in *sandbox*

To do:  

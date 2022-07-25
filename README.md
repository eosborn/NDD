# Variant and Gene Identification in Neurodevelopmental Delay Cohort 
## project for the Translational Health Data Science Fellowship

### Aim 1: collect list of genes associated with neurodevelopemntal disorders (ASD, speech delay, nonverbal)
- see `NDD_Gene_List.txt` in */Aim_1* directory
- generated using a broad literature search

To do:
	1. add literature sources to repo 
	
### Aim 2: scrape KEGG mapper output for pathways/modules and create dictionaires to organize groupings
- raw KEGG mapper output located in *Aim_2/* as .txt files
- parsed output (just grouping names and gene names) located in *Aim_2/parsed/*
- dictionaries with gene names as *keys* located in *Aim_2/dicts/gene_dicts/*
	- made using `dict_gene_script.py` python script
- dictionaries with grouping names as *keys* located in *Aim_2/dicts/grouping_dicts* 
	- made using `dict_grouping_script.py` python script
	
To do:
	1. include shell script used to parse KEGG output prior to dictionary creation
	2. add *README* to directory

### Aim 3: apply gene list and groupings to existing sequencing data for patient cohort
- see `add_groups.py` and `counting.py` python scripts located in *Aim_3/* directory
To do:
	1. add *README* to directory 
	2. fill in function docstrings in `counting.py`

### Aim 4: identify relevant genes/groupings for clinical interpretation  
- currently workshopping in *sandbox*

To do:
	1. 

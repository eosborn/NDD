# NDD
## project for the Translational Health Data Science Fellowship

###Aim 1: collect list of genes associated with neurodevelopemntal disorders (ASD, speech delay, nonverbal)
- see `NDD_Gene_List.txt`

###Aim 2: scrape KEGG mapper output for pathways/modules and create dictionaires to organize groupings
- raw output located in *KEGG_output/*
- parsed output (just grouping names and gene names) located in *KEGG_output/parsed*
- dictionaries with gene names as *keys* located in *KEGG_output/dicts/gene_dicts*
	- made using `dict_gene_script.py` python script
- dictionaries with grouping names as *keys* located in *KEGG_output/dicts/grouping_dicts* 
	- made using `dict_grouping_script.py` python script

###Aim 3: apply gene list and groupings to existing sequencing data for patient cohort, identify relevant genes/groupings for clinical interpretation  

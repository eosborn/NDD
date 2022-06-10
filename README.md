# Variant and Gene Identification in Neurodevelopmental Delay Cohort 
## project for the Translational Health Data Science Fellowship

### Aim 1: collect list of genes associated with neurodevelopemntal disorders (ASD, speech delay, nonverbal)
- see `NDD_Gene_List.txt`
- generated using a broad literature search

To do:
	1. add literature sources to repo 
	
### Aim 2: scrape KEGG mapper output for pathways/modules and create dictionaires to organize groupings
- raw output located in *KEGG_output/*
- parsed output (just grouping names and gene names) located in *KEGG_output/parsed*
- dictionaries with gene names as *keys* located in *KEGG_output/dicts/gene_dicts*
	- made using `dict_gene_script.py` python script
- dictionaries with grouping names as *keys* located in *KEGG_output/dicts/grouping_dicts* 
	- made using `dict_grouping_script.py` python script

### Aim 3: apply gene list and groupings to existing sequencing data for patient cohort
- see `add_groups.py` python script located in *Filtering/*

To do:
	1. make `add_groups.py` more efficient --> input list of files to loop over so that dictionaries don't need to be recreated every time the script is called 
	2. make `add_groups.py` executable from the command line

### Aim 4: identify relevant genes/groupings for clinical interpretation  

To do:
	1. due to discovery of artifacts in original patient VCF files used for filtering, re-align raw FASTQ files, re-create VCF files and re-run filtering pipeline (including steps outlined in Aims 1-3)
	2. write a script that groups samples by variants occuring in common groupings (pathways, networks, etc) 
	3. using IGV and mapped+annotated BAM files, confirm that filtered variants of interest are real

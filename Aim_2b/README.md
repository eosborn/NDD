
# Aim 2b

### add_groups.py
goal: integrate KEGG groupings into cohort vcf files

to use:
`python add_goruops.py <vcf_file_directory> <grouping_1_name> <grouping_1_dict_file> <grouping_2_name> ... `

approach:
1. for each group, create a new column in the vcf files
2. for each row of the vcf file (each row = one variant), look up that variant's gene's associated groups in grouping dictionary and add to vcf file
3. write out appended vcf files as both txt and feather files

### counting_v2.py
goal: create files that count genes and groupings across all input feather-formatted vcf files

to use:
`python counting_v2.py <vcf_file_directory> <out_directory> <meta_data_directory>` 

approach:
1. create empty dataframes to input counting data
2. for each vcf file, count the number of occurences of genes and groupings, and input counts into respective count dfs
3. append metadata to count dfs 
4. write out count dfs as txt and pickle files

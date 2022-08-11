#!/bin/bash

read -p 'Enter files: ' files
#echo $files

mkdir parsed_to_tsv

#create genes_tsv.txt
mkdir parsed_to_tsv/col_files
mkdir parsed_to_tsv/grouping_genes_files 
for file in $files
do
	echo $file
	read -p 'Enter grouping name: ' name
	#create seperate files ("Argument list too long" error if tried to save to variables)
	awk '$1=$1' $file | grep "hsa:" | cut -f2 -d' ' | \
		sed 's/;//' > parsed_to_tsv/col_files/"${name}_col1.txt"  
	awk '$1=$1' $file | grep "hsa:" | cut -f1 -d' ' | \
		sed 's/;//' > parsed_to_tsv/col_files/"${name}_col2.txt" 
	awk '$1=$1' $file | grep "hsa:" | sed 's/^[^;]*;//' | \
                sed 's/[[:space:]]//' >  parsed_to_tsv/col_files/"${name}_col3.txt" 
	paste parsed_to_tsv/col_files/"${name}_col1.txt" parsed_to_tsv/col_files/"${name}_col2.txt" \
		parsed_to_tsv/col_files/"${name}_col3.txt" | sort | uniq > \
		parsed_to_tsv/grouping_genes_files/"${name}_genes_tsv.txt"

	#create grouping tsv
	cat $file | grep -v 'hsa:' | awk '$1=$1' | cut -f1 -d' ' \
		> parsed_to_tsv/col_files/"only_${name}_col1.txt"
	cat $file | grep -v 'hsa:' | awk '$1=$1' | sed 's/[^ ]*[[:space:]]//' | \
		sed 's/- Homo sapiens (human)[[:space:]]//' | sed 's/[[:space:]](.*//' \
		> parsed_to_tsv/col_files/"only_${name}_col2.txt"
	cat $file | grep -v 'hsa:' | awk '$1=$1' | sed 's/.*(//' | sed 's/)//' \
		> parsed_to_tsv/col_files/"only_${name}_col3.txt"
	paste parsed_to_tsv/col_files/"only_${name}_col1.txt"  \
		parsed_to_tsv/col_files/"only_${name}_col2.txt" \
		parsed_to_tsv/col_files/"only_${name}_col3.txt" \
		> parsed_to_tsv/"only_${name}_tsv.txt"

done
cat parsed_to_tsv/grouping_genes_files/*txt | sort | uniq > parsed_to_tsv/genes_tsv.txt

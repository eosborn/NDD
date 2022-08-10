#!/bin/bash

read -p 'Enter files: ' files
#echo $files

mkdir parsed_to_tsv
mkdir parsed_to_tsv/col_files 

for file in $files
do
	echo $file
	read -p 'Enter grouping name: ' name
	#create seperate files ("Argument list too long" error if tried to save to variables)
	awk '$1=$1' $file | grep "hsa:" | cut -f2 -d' ' | \
		sed 's/;//' | uniq > parsed_to_tsv/col_files/"${name}_col1.txt"  
	awk '$1=$1' $file | grep "hsa:" | cut -f1 -d' ' | \
		sed 's/;//' | uniq > parsed_to_tsv/col_files/"${name}_col2.txt" 
	awk '$1=$1' $file | grep "hsa:" | sed 's/^[^;]*;//' | \
                sed 's/[[:space:]]//' | uniq > parsed_to_tsv/col_files/"${name}_col3.txt" 
	paste parsed_to_tsv/col_files/"${name}_col1.txt" parsed_to_tsv/col_files/"${name}_col2.txt" \
		parsed_to_tsv/col_files/"${name}_col3.txt" | sort | uniq > parsed_to_tsv/"${name}_tsv.txt"

done

cat parsed_to_tsv/*txt | sort | uniq > parsed_to_tsv/genes_tsv.txt

#!/bin/bash

#input grouped.txt variant files
read -p 'Enter files: ' files
#files=test_data2/*grouped.txt
path=`echo $files[1] | sed s/"\/[^ ]*"//`

read -p 'Enter gene list: ' genes
#genes=../Aim_3/trunc_mult.txt

mkdir trunc_vars

for file in $files
do 
	sample=`echo $file | sed s/$path"\/"// | cut -f1 -d'_'`
	echo $sample >> trunc_vars/all_trunc_vars.txt
	
	cat $file | grep -f $genes | grep -vw 'HTT' | \
		grep -vw 'nonframeshift\|ExonicFunc\|nonsynonymous' | \
		tee -a trunc_vars/all_trunc_vars.txt trunc_vars/${sample}_trunc_vars.txt
done


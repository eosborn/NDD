### Aim 3

#### *all_sub_genes.txt*
all genes appearing in cohort, excluding *TTN*, *MUC12*, and *FLG* which were flagged to be not be of interest due to mismapping
- 1102 genes
> cut -f1 -d',' ~/Aim_2b/cohort_count_output/genes_subset.csv | grep -v 'gene' > all_sub_genes.txt

#### *trunc_all.txt*
all genes (excluding flagged genes listed above) containing a truncating mutation 
- 140 genes
- *grouped.txt variant files not currently located in repo*
> cat *grouped.txt | grep -vw 'TTN\|MUC12\|FLG' | grep -vw 'nonframeshift\|ExonicFunc\|nonsynonymous' | cut -f7 | sort | uniq
 > trunc_all.txt

#### *trunc_mult.txt*
genes (excluding flagged genes) containing a multiple truncating mutations across the cohort 
- 31 genes
> cat *grouped.txt | grep -vw 'TTN\|MUC12\|FLG' | grep -vw 'nonframeshift\|ExonicFunc\|nonsynonymous' | cut -f7 | sort | uniq -c | grep -vw '1' | sed 's/.*[[:space:]]//' > trunc_mult.txt

#### *trunc_mult_counts.txt*
trunc_mult genes with second count column (# of times that gene has a truncating mutation across the cohort)
> cat *grouped.txt | grep -vw 'TTN\|MUC12\|FLG' | grep -vw 'nonframeshift\|ExonicFunc\|nonsynonymous' | cut -f7 | sort | uniq -c | grep -vw '1' | awk '$1=$1' | cut -f1 -d' ' > tmp.txt
> paste trunc_mult.txt tmp.txt > trunc_mult_counts.txt
> rm tmp.txt

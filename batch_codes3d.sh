#!/bin/bash

#snp_file="/mnt/3dgenome/projects/tfad334/hvn/diabetes/results/codes3d_results/snps/undone_snps.txt"
#wc -l $snp_file
snpdir="../GWAS-Catalog-Analysis/data/batch_snps/done"
#split -l 50 -d $snp_file
#mv $snp_file ../
#mkdir ../data/batch_results1
#mkdir ../data/batch_snps/done1
results_dir="../GWAS-Catalog-Analysis/data/batch_results1"
done_dir="../GWAS-Catalog-Analysis/data/batch_snps/done1"
for f in $snpdir/*
do
    echo Processing $f 
    python codes3d/codes3d.py -i $f -o $results_dir/$f -r mboi 
    mv $f done_dir
done

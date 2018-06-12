#!/bin/bash

#snp_file="/mnt/3dgenome/projects/tfad334/all-gwas/data/snps/all_catalog_snps.txt"
#wc -l $snp_file

#split -l 50 $snp_file
#mv $snp_file ../
#mkdir /mnt/3dgenome/projects/tfad334/all-gwas/data/snps/results
for f in /mnt/3dgenome/projects/tfad334/all-gwas/data/snps/*
do
    echo Processing $f
    python ./codes3d.py -i $f -o mnt/3dgenome/projects/tfad334/all-gwas/data/snps/results/$f -x GM12878_TEST -p 4 -c ../docs/conf.py

done
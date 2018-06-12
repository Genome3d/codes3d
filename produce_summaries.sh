#!/bin/bash


#repdir=/mnt/3dgenome/projects/tfad334/hvn/analysis/control/random_daSNPs
repdir=/mnt/3dgenome/projects/tfad334/GWAS-Catalog-Analysis/results/undone
#dirs= ls -dtr $repdir/*/ | less -500

#for i in $(ls -dtr $repdir/*/);
for i in $(ls $repdir/);
do
    echo 'Producing summary for...'
    #echo -e '\t ...'$i
    #python codes3d/produce_summary.py -e $i/eqtls.txt -o $i  -c docs/codes3d.conf
	python codes3d/produce_summary.py -e $repdir/$i/eqtls.txt  -o $repdir/$i/ -c docs/codes3d.conf
	#if [ ! -f `$repdir/$i/significant_eqtls.txt`] 
	#then
	#	echo "Done!"
	#fi
done


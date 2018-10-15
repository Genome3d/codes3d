#!/usr/bin/env python

import argparse,codes3d,configparser,os


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Helper script for building or updating the SNP index from a directory containing SNP files in BED format.")
	parser.add_argument("-s","--snp_dir",required=True,help="The directory containing the SNP files in BED format.")
	parser.add_argument("-i","--id_col",help="The column containing the ID of the SNP (default is column 4 in dbSNP files).")
	parser.add_argument("-c","--chr_col",help="The column containing the chromosome of the SNP (default is column 1 in dbSNP files).")
	parser.add_argument("-l","--locus_col",help="The column containing the locus of the SNP on the chromosome (\"start\" is used by default, column 2 in dbSNP files).")
	parser.add_argument("-C","--config_file",default=os.path.join(os.path.dirname(__file__),"../docs/codes3d.conf"),help="The configuration file to be used in this instance (default: conf.py)")
	parser.add_argument("-o","--output_fp",help="The file to which to output the the results (default: whatever SNP_DATABASE_FP is specified in the config file supplied).")
	args = parser.parse_args()

	config = configparser.ConfigParser()
	config.read(args.config_file)
	codes3d.build_snp_index(args.snp_dir,args.output_fp,config,args.id_col,args.chr_col,args.locus_col)
	

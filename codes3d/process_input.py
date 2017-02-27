#!/usr/bin/env python

from configobj import ConfigObj
import argparse,codes3d,os

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',required=True,help="The the dbSNP IDs or loci of SNPs of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-c","--config",default="docs/conf.py",help="The configuration file to be used in this instance (default: conf.py)")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	args = parser.parse_args()
	config = ConfigObj(args.config)
	snp_database_fp = config["SNP_DATABASE_FP"]
	fragment_database_fp = config["FRAGMENT_DATABASE_FP"]
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)

	codes3d.process_inputs(args.inputs,snp_database_fp,fragment_database_fp,args.output_dir)
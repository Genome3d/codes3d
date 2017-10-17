#!/usr/bin/env python

import argparse,codes3d,configparser,os

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-e","--eqtl_file",help="The list of eQTLs generated by find_eqtls.py")
	parser.add_argument("-c","--config",default="docs/conf.py",help="The configuration file to be used in this instance (default: conf.py)")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	args = parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config)
	expression_table_fp = config.get("Defaults","EXPRESSION_TABLE_FP")
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)
	
	eqtls = codes3d.parse_eqtls_file(args.eqtl_file)
	gene_exp = codes3d.get_gene_expression_information(eqtls,expression_table_fp,args.output_dir)
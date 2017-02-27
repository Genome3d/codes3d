#!/usr/bin/env python

import argparse,codes3d


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Script for indexing eQTL tables (defaults are for tables of significant eQTLs from GTEx). NB: ALL INDICES ARE 1-BASED, NOT 0-BASED.")
	parser.add_argument("-e","--input_eqtl_table_fp",required=True,help="The eQTL table to be indexed.")
	parser.add_argument("-s","--snp_id_column",help="The column containing the SNP identifier (default: 23, rsID column in GTEx tables.)")
	parser.add_argument("-g","--gene_symbol_column",help="The column containing the gene symbol (e.g. \"CELSR2\", default: 27.)")
	parser.add_argument("-c","--gene_chr_column",help="The column containing the gene chromosome (default: 29.)")
	parser.add_argument("-b","--gene_start_column",help="The column containing the start coordinate of the gene (default: 30.)")
	parser.add_argument("-f","--gene_stop_column",help="The column containing the end coordinate of the gene (default: 31.)")
	parser.add_argument("-p","--p_val_column",help="The column containing the p-value for the eQTL (deafult: 6.)")
	parser.add_argument("-x","--effect_size_column",help="The column containing the effect size of the eQTL (default: 3.)")
	parser.add_argument("-o","--output_fp",help="The output file for the new eQTL index (default: same as the filename of the input file with extension \".db\")")
	args = parser.parse_args()
	
	codes3d.build_eqtl_table_index(args.input_eqtl_table_fp,args.snp_id_column,args.gene_symbol_column,args.gene_chr_column,args.gene_start_column,args.gene_stop_column,args.p_val_column,args.effect_size_column,args.output_fp)
	

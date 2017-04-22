#!/usr/bin/python

import argparse,codes3d,configparser


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Create a BED file detailing the locations of genes in the genome, and a database containing additional gene information. Note: If a file in .gtf format is supplied, no other arguments are required.")
	parser.add_argument("-i","--gene_files",required=True,nargs='+',help="The gene file/s to be indexed; either in tabular format, or, by default, the .gtf file format, as supplied by the GTEx project.")
	parser.add_argument("-g","--symbol_col",type=int,help="The index of the column containing the gene symbol (non-zero based; default: ).")
	parser.add_argument("-c","--chr_col",type=int,help="The index of the column containing the chromosome name (non-zero based; default: ).")
	parser.add_argument("-s","--start_col",type=int,help="The index of the column containing the gene start site (non-zero based; default: ).")
	parser.add_argument("-e","--end_col",type=int,help="The index of the column containing the gene end site (non-zero based; default: ).")
	parser.add_argument("-p","--p_threshold_col",type=int,help="The index of the column containing the GTEx p-threshold for this gene (optional; non-zero based; default: ).")
	parser.add_argument("-H","--no_header",action="store_true",help="Use this option if the table has no header.")
	parser.add_argument("-b","--output_bed_fp",help="The path to which to output the resultant BED file of gene locations (default: the input file name with the extension \".bed\").")
	parser.add_argument("-o","--output_db",help="The path to which to output the resultant gene index database (default: the input file name with the extension \".db\").")
	parser.add_argument("-C","--config_file",default=os.path.join(os.path.dirname(__file__),"../docs/codes3d.conf"),help="The configuration file specifying the location of the CoDeS3D library (default: docs/codes3d.conf).")
	args = parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config_file)

	build_gene_index(args.gene_files,args.output_bed_fp,args.output_db,config,args.symbol_col,args.chr_col,args.start_col,args.end_col,args.p_threshold_col,args.no_header)
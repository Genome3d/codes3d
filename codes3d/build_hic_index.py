#!/usr/bin/env python

import argparse,codes3d


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Tool for indexing HiC tables into databases for use by the CoDeS3D toolkit.")
	parser.add_argument("-i","--input_hic_fp",required=True,help="The filepath of the table to index.")
	parser.add_argument("-c1","--chr1_col",type=int,help="Column number of partner 1 chromosome (default: 3, NOT ZERO-BASED)")
	parser.add_argument("-c2","--chr2_col",type=int,help="Column number of partner 2 chromosome (default: 5, NOT ZERO-BASED)")
	parser.add_argument("-f1","--frag1_col",type=int,help="Column number of partner 1 fragment (default: 7, NOT ZERO-BASED)")
	parser.add_argument("-f2","--frag2_col",type=int,help="Column number of partner 2 fragment (default: 9, NOT ZERO-BASED)")
	parser.add_argument("-m1","--mapq1_col",type=int,help="Column number of partner 1 MAPQ score (default: 10, NOT ZERO-BASED)")
	parser.add_argument("-m2","--mapq2_col",type=int,help="Column number of partner 2 MAPQ score (default: 11, NOT ZERO-BASED)")
	parser.add_argument("-m","--mapq_cutoff",type=int,help="The minimum MAPQ score for a read to be considered (default: 150)")
	parser.add_argument("-o","--output_fp",help="The output file for the new HiC index (default: same as the filename of the input file with extension \".db\")")
	args = parser.parse_args()

	codes3d.build_hic_index(args.input_hic_fp,args.output_fp,args.chr1_col,args.chr2_col,args.frag1_col,args.frag2_col,args.mapq1_col,args.mapq2_col,args.mapq_cutoff)

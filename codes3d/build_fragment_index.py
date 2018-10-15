#!/usr/bin/env python

import argparse,codes3d


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Tool for indexing fragments from fragmented genome file.")
	parser.add_argument("-f","--fragment_fp",required=True,help="The fragmented genome file to index.")
	parser.add_argument("-o","--output_fp",help="The file to which to output the resultant index and BED file (default: the name of the input file, with the extension \".db\"/\".bed\", respectively.)")
	args = parser.parse_args()

	codes3d.build_fragment_index(args.fragment_fp,args.output_fp)
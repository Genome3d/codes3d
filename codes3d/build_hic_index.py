from sets import Set
import argparse,sqlite3,os

def build_hic_table_index(input_hic_fp,chr1_col,chr2_col,frag1_col,frag2_col,mapq1_col,mapq2_col,mapq_cutoff,output_fp):
        ##Do line count for progress meter
	print "Determining table size..."
	hic_table = open(input_hic_fp,'r')
	lines = 0
	for i in hic_table:
	    lines += 1
	hic_table.close()
	lines = lines//100*100 #Get an approximation
	do_linecount = not lines == 0
	
	int_db = sqlite3.connect(output_fp)
	interactions = int_db.cursor()
	interactions.execute("CREATE TABLE IF NOT EXISTS interactions (chr1 text, fragment1 text, chr2 text, fragment2 text)")
	interactions.execute("CREATE INDEX IF NOT EXISTS i_1 ON interactions (chr1, fragment1)")
	chr1_col = chr1_col - 1
	chr2_col = chr2_col - 1
	frag1_col = frag1_col - 1
	frag2_col = frag2_col - 1
	mapq1_col = mapq1_col - 1
	mapq2_col = mapq2_col - 1
	with open(input_hic_fp,'r') as rao_table:
		print "Indexing HiC interaction table..."
		for i,line in enumerate(rao_table):
			if do_linecount:
				if i % (lines/100) == 0:
					print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
			interaction = line.strip().split(' ')
			if int(interaction[mapq1_col]) >= mapq_cutoff and int(interaction[mapq2_col]) >= mapq_cutoff:
				chr1 = interaction[chr1_col]
				frag1 = interaction[frag1_col]
				chr2 = interaction[chr2_col]
				frag2 = interaction[frag2_col]
				interactions.execute("INSERT INTO interactions VALUES(?,?,?,?)", [chr1,frag1,chr2,frag2])
				interactions.execute("INSERT INTO interactions VALUES(?,?,?,?)", [chr2,frag2,chr1,frag1])
        int_db.commit()
        interactions.close()
	print "Done indexing HiC interaction table."

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Tool for indexing HiC tables into databases for use by the CoDeS3D toolkit.")
	parser.add_argument("-i","--input_hic_fp",required=True,help="The filepath of the table to index.")
	parser.add_argument("-c1","--chr1_col",type=int,default=3,help="Column number of partner 1 chromosome (default: 3, NOT ZERO-BASED)")
	parser.add_argument("-c2","--chr2_col",type=int,default=7,help="Column number of partner 2 chromosome (default: 5, NOT ZERO-BASED)")
	parser.add_argument("-f1","--frag1_col",type=int,default=5,help="Column number of partner 1 fragment (default: 7, NOT ZERO-BASED)")
	parser.add_argument("-f2","--frag2_col",type=int,default=9,help="Column number of partner 2 fragment (default: 9, NOT ZERO-BASED)")
	parser.add_argument("-m1","--mapq1_col",type=int,default=10,help="Column number of partner 1 MAPQ score (default: 10, NOT ZERO-BASED)")
	parser.add_argument("-m2","--mapq2_col",type=int,default=11,help="Column number of partner 2 MAPQ score (default: 11, NOT ZERO-BASED)")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="The minimum MAPQ score for a read to be considered (default: 150)")
	parser.add_argument("-o","--output_fp",help="The output file for the new HiC index (default: same as the filename of the input file with extension \".db\")")
	args = parser.parse_args()
	if not args.output_fp:
		if not args.input_hic_fp.rfind('.') == -1:
			args.output_fp = args.input_hic_fp[:args.input_hic_fp.rfind('.')] + ".db"
		else:
			args.output_fp = args.input_hic_fp + ".db"

	build_hic_table_index(args.input_hic_fp,args.chr1_col,args.chr2_col,args.frag1_col,args.frag2_col,args.mapq1_col,args.mapq2_col,args.mapq_cutoff,args.output_fp)

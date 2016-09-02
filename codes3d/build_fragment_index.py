import argparse,sqlite3,os

def build_fragment_index(fragment_fp,output_fp):
	if os.path.isfile(output_fp):
		os.remove(output_fp)
	fragment_index_db = sqlite3.connect(output_fp)
	fragment_index = fragment_index_db.cursor()
	fragment_index.execute("CREATE TABLE fragments (chr text, start integer, end integer, fragment integer)")
	fragment_index.execute("CREATE INDEX f_index ON fragments (chr,fragment)")
	
	with open(fragment_fp,'r') as fragments_bed:
		for line in fragments_bed:
			fragment = line.strip().split('\t')
			fragment_index.execute("INSERT INTO fragments VALUES (?,?,?,?)", [fragment[0][fragment[0].find("chr")+3:],int(fragment[1]),int(fragment[2]),fragment[3]])
	fragment_index_db.commit()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Tool for indexing fragments from fragmented genome file.")
	parser.add_argument("-f","--fragment_fp",required=True,help="The fragmented genome file to index.")
	parser.add_argument("-o","--output_fp",help="The file to which to output the resultant index (default: the name of the input file, with the extension \".db\")")
	args = parser.parse_args()
	if not args.output_fp:
		if not args.fragment_fp.rfind('.') == -1:
			args.output_fp = args.fragment_fp[:args.fragment_fp.rfind('.') == -1] + ".db"
		else:
			args.output_fp = args.fragment_fp + ".db"

	build_fragment_index(args.fragment_fp,args.output_fp)
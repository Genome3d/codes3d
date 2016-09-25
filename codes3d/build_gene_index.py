#!/usr/bin/python
import argparse,sqlite3,os

def build_gene_index(gene_files,symbol_col,chr_col,start_col,end_col,p_thresh_col,no_header,output_bed,output_db):
	genes = {}
	symbol_col -= 1
	chr_col -= 1
	start_col -= 1
	end_col -= 1
	if p_thresh_col:
		p_thresh_col -= 1
	if os.path.isfile(output_db):
		upsert = raw_input("WARNING: Upserting input to existing SNP database (%s). Continue? [y/N] " % output_db)
		if not upsert.lower() == 'y':
			print "Did not write to existing SNP database."
			return
	gene_index_db = sqlite3.connect(output_db)
	gene_index = gene_index_db.cursor()
	if p_thresh_col:
		gene_index.execute("CREATE TABLE IF NOT EXISTS genes (symbol text, chr text, start integer, end integer, double p_thresh)")
	else:
		gene_index.execute("CREATE TABLE IF NOT EXISTS genes (symbol text, chr text, start integer, end integer)")
	gene_index.execute("CREATE INDEX IF NOT EXISTS g_index ON genes (symbol)")
	for gene_file in gene_files:
		with open(gene_file,'r') as genefile:
			if not no_header:
				genefile.readline()
			line = genefile.readline()
			while line:
				gene = line.strip().split('\t')
				gene_symbol = gene[symbol_col]
				if not genes.has_key(gene_symbol):
					genes[gene_symbol] = { "start": int(gene[start_col]), "end": int(gene[end_col]) }
					if gene[chr_col].startswith("chr"):
						genes[gene_symbol]["chr"] = gene[chr_col][3:]
					else:
						genes[gene_symbol]["chr"] = gene[chr_col]
					if p_thresh_col:
						try:
							genes[gene_symbol]["p_thresh"] = float(gene[p_thresh_col])
						except ValueError:
							genes[gene_symbol]["p_thresh"] = gene[p_thresh_col]
				else:
					curr_length = genes[gene_symbol]["end"] - genes[gene_symbol]["start"]
					start = int(gene[start_col])
					end = int(gene[end_col])
					if curr_length < (end - start):
						genes[gene_symbol]["start"] = start
						genes[gene_symbol]["end"] = end
				line = genefile.readline()
	bed_out = None
	if output_bed:
		bed_out = open(output_bed,'w')
	for gene in genes.keys():
		if bed_out:
			bed_out.write('chr%s\t%s\t%s\t%s\n' % (genes[gene]["chr"],genes[gene]["start"],genes[gene]["end"],gene))
		if p_thresh_col:
			gene_index.execute("INSERT INTO genes VALUES (?,?,?,?,?)", [gene,genes[gene]["chr"],genes[gene]["start"],genes[gene]["end"],genes[gene]["p_thresh"]])
		else:
			gene_index.execute("INSERT INTO genes VALUES (?,?,?,?)", [gene,genes[gene]["chr"],genes[gene]["start"],genes[gene]["end"]])
	gene_index_db.commit()
	gene_index.close()
	if bed_out:
		bed_out.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--gene_files",required=True,nargs='+')
	parser.add_argument("-g","--symbol_col",type=int,default=27)
	parser.add_argument("-c","--chr_col",type=int,default=29)
	parser.add_argument("-s","--start_col",type=int,default=30)
	parser.add_argument("-e","--end_col",type=int,default=31)
	parser.add_argument("-p","--p_threshold_col",type=int)
	parser.add_argument("-H","--no_header",action="store_true",default=False)
	parser.add_argument("-b","--output_bed_fp")
	parser.add_argument("-o","--output_db")
	args = parser.parse_args()
	if not args.output_db:
		args.output_db = "gene_index.db"
	
	build_gene_index(args.gene_files,args.symbol_col,args.chr_col,args.start_col,args.end_col,args.p_threshold_col,args.no_header,args.output_bed_fp,args.output_db)
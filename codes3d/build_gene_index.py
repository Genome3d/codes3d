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
		#Do line count for progress meter
		lines = 0
		with open(gene_file,'r') as genefile:
			print "Determining table size..."
			for i in genefile:
				lines += 1
			lines = lines//100*100 #Get an approximation
			do_linecount = not lines == 0
		
		with open(gene_file,'r') as genefile:
			#Determine if input file is a GTEx-supplied gene reference
			is_gtex_file = gene_file.endswith(".gtf")
			#If so, process headers accordingly
			if is_gtex_file:
				line = genefile.readline()
				while line.startswith("##"):
					line = genefile.readline()
			#Otherwise, process headers according to script options
			else:
				if not no_header:
					genefile.readline()
				line = genefile.readline()
			i = 0
			#For each line
			while line:
				if do_linecount:
					if i % (lines/100) == 0:
						print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
				#If the line is a GTEx file, extract information accordingly
				if is_gtex_file:
					entry = line.strip().split('\t')
					if entry[2] == "gene":
						gene_stats = entry[8].split(';')
						gene_symbol = gene_stats[4].strip().split(' ')[1].strip('"')
						gene_chr = entry[0]
						gene_start = int(entry[3])
						gene_end = int(entry[4])
					else:
						i += 1
						line = genefile.readline()
						continue #Skip if the entry is not for a canonical gene
				#Otherwise, extract by column
				else:			
					gene = line.strip().split('\t')
					gene_symbol = gene[symbol_col]
					if gene[chr_col].startswith("chr"):
						gene_chr = gene[chr_col][3:]
					else:
						gene_chr = gene[chr_col]
					gene_start = int(gene[start_col])
					gene_end = int(gene[end_col])
					if p_thresh_col:
						try:
							gene_p_thresh = float(gene[p_thresh_col])
						except ValueError:
							gene_p_thresh = gene[p_thresh_col]
				#Enter into index, regardless of input file type
				if not genes.has_key(gene_symbol):
					genes[gene_symbol] = { "chr": gene_chr, "start": gene_start, "end": gene_end }
					if p_thresh_col:
						genes[gene_symbol]["p_thresh"] = gene_p_thresh
				else:
					curr_length = genes[gene_symbol]["end"] - genes[gene_symbol]["start"]
					if curr_length < abs(gene_end - gene_start):
						genes[gene_symbol]["start"] = gene_start
						genes[gene_symbol]["end"] = gene_end
				line = genefile.readline()
				i += 1
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
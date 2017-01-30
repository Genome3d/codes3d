import argparse,os,sqlite3

def build_eqtl_table_index(table_fp,snp_col,gene_symbol_col,gene_chr_col,gene_start_col,gene_stop_col,p_val_col,effect_size_col,output_fp):
	if os.path.isfile(output_fp):
		upsert = raw_input("WARNING: Upserting input to existing eQTL database %s. Continue? [y/N] " % output_fp)
		if not upsert.lower() == 'y':
			print "Did not write to existing eQTL database."
			return

	snp_col -= 1
	gene_symbol_col -= 1
	gene_chr_col -= 1
	gene_start_col -= 1
	gene_stop_col -= 1

	effect_size_col -= 1
	table_index_db = sqlite3.connect(output_fp)
	table_index = table_index_db.cursor()
	table_index.execute("CREATE TABLE IF NOT EXISTS eqtls (rsID text, gene_name text, gene_chr text, gene_start integer, gene_end integer, pvalue real, effect_size real)")
	table_index.execute("CREATE INDEX IF NOT EXISTS id ON eqtls (rsID)")
	##Do line count for progress meter
	do_linecount = True
	print "Determining table size..."
	with open(table_fp,'r') as eqtl_table
		lines = 0
		for i in eqtl_table:
			lines += 1
	lines = lines//100*100 #Get an approximation
	do_linecount = not lines == 0
	
	with open(table_fp,'r') as eqtl_table:
		for i,line in enumerate(eqtl_table):
			if do_linecount:
				if i % (lines/100) == 0:
					print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
			if i == 0:
				continue
			eqtl = line.strip().split('\t')
			table_index.execute("INSERT INTO eqtls VALUES (?,?,?,?,?,?,?)",[eqtl[snp_col],eqtl[gene_symbol_col],eqtl[gene_chr_col],eqtl[gene_start_col],eqtl[gene_stop_col],eqtl[p_val_col],eqtl[effect_size_col]])
	table_index_db.commit()
	print "Done indexing eQTL table."

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Script for indexing eQTL tables (defaults are for tables of significant eQTLs from GTEx). NB: ALL INDICES ARE 1-BASED, NOT 0-BASED.")
	parser.add_argument("-e","--input_eqtl_table_fp",required=True,help="The eQTL table to be indexed.")
	parser.add_argument("-s","--snp_id_column",default=23,help="The column containing the SNP identifier (default: 23, rsID column in GTEx tables.)")
	parser.add_argument("-g","--gene_symbol_column",default=27,help="The column containing the gene symbol (e.g. \"CELSR2\", default: 27.)")
	parser.add_argument("-c","--gene_chr_column",default=29,help="The column containing the gene chromosome (default: 29.)")
	parser.add_argument("-b","--gene_start_column",default=30,help="The column containing the start coordinate of the gene (default: 30.)")
	parser.add_argument("-f","--gene_stop_column",default=31,help="The column containing the end coordinate of the gene (default: 31.)")
	parser.add_argument("-p","--p_val_column",default=6,help="The column containing the p-value for the eQTL (deafult: 6.)")
	parser.add_argument("-x","--effect_size_column",default=3,help="The column containing the effect size of the eQTL (default: 3.)")
	parser.add_argument("-o","--output_fp",help="The output file for the new eQTL index (default: same as the filename of the input file with extension \".db\")")
	args = parser.parse_args()
	if not args.output_fp:
		if not args.input_eqtl_table_fp.rfind('.') == -1:
			args.output_fp = args.input_eqtl_table_fp[:args.input_eqtl_table_fp.rfind('.')] + ".db"
		else:
			args.output_fp = args.input_eqtl_table_fp + ".db"

	build_eqtl_table_index(args.input_eqtl_table_fp,args.snp_id_column,args.gene_symbol_column,args.gene_chr_column,args.gene_start_column,args.gene_stop_column,args.p_val_column,args.effect_size_column,args.output_fp)
	
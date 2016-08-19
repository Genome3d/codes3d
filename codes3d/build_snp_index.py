from configobj import ConfigObj
import argparse,sqlite3,os

def build_snp_index(snp_dir,id_col,chr_col,locus_col,output_fp):
	print "Building SNP index..."
	if not os.path.isdir(snp_dir):
		print "Error: argument to build SNP index must be a directory."
		return
	if os.path.isfile(output_fp):
		upsert = raw_input("WARNING: Upserting input to existing SNP database (%s). Continue? [y/N] " % output_fp)
		if not upsert.lower() == 'y':
			print "Did not write to existing SNP database."
			return
	snp_index_db = sqlite3.connect(output_fp)
	snp_index = snp_index_db.cursor()
	snp_index.execute("CREATE TABLE IF NOT EXISTS snps (rsID text unique, chr text, locus integer)")
	snp_index.execute("CREATE INDEX IF NOT EXISTS id ON snps (rsID,chr,locus)")
	res = ""
	id_col = id_col - 1
	chr_col = chr_col - 1
	locus_col = locus_col - 1
	for bed_fp in os.listdir(snp_dir):
		print "\tProcessing " + bed_fp
		bed = open(snp_dir + '/' + bed_fp,'r')
		for line in bed:
			if not line.startswith("chr"):
				continue
			snp = line.strip().split('\t')
			try:
				snp_index.execute("INSERT INTO snps VALUES(?,?,?)",[snp[id_col],snp[chr_col][snp[chr_col].find("chr")+3:],int(snp[locus_col])])
			except sqlite3.IntegrityError:
				if res == "":
					res = raw_input("Warning: database already contains an entry for SNP with ID %s. Overwrite?\n1: Yes\n2: No (default)\n3: Yes To All\n4: No To All\nChoice: " % snp[id_col])
				if res.strip() == "1":
					print "Overwriting SNP %s" % snp[id_col]
					snp_index.execute("DELETE FROM snps WHERE rsID=?",[snp[id_col]])
					snp_index.execute("INSERT INTO snps VALUES(?,?,?)",[snp[id_col],snp[chr_col][snp[chr_col].find("chr")+3:],int(snp[locus_col])])
					res = ""
				elif res.strip() == "3":
					print "Overwriting SNP %s" % snp[id_col]
					snp_index.execute("DELETE FROM snps WHERE rsID=?",[snp[id_col]])
					snp_index.execute("INSERT INTO snps VALUES(?,?,?)",[snp[id_col],snp[chr_col][snp[chr_col].find("chr")+3:],int(snp[locus_col])])
				elif res.strip() == "4":
					print "Skipping input SNP %s" % snp[id_col]
					pass
				else:
					print "Skipping input SNP %s" % snp[id_col]
					res = ""
		bed.close()
	print "\tWriting SNP index to file..."
	snp_index_db.commit()
	print "Done building SNP index."

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Helper script for building or updating the SNP index from a directory containing SNP files in BED format.")
	parser.add_argument("-s","--snp_dir",required=True,help="The directory containing the SNP files in BED format.")
	parser.add_argument("-i","--id_col",default=4,help="The column containing the ID of the SNP (default is column 4 in dbSNP files).")
	parser.add_argument("-c","--chr_col",default=1,help="The column containing the chromosome of the SNP (default is column 1 in dbSNP files).")
	parser.add_argument("-l","--locus_col",default=2,help="The column containing the locus of the SNP on the chromosome (\"start\" is used by default, column 2 in dbSNP files).")
	parser.add_argument("-C","--config",default="docs/conf.py",help="The configuration file to be used in this instance (default: conf.py)")
	parser.add_argument("-o","--output_fp",help="The file to which to output the the results (default: whatever SNP_DATABASE_FP is specified in the config file supplied).")
	args = parser.parse_args()
	if not args.output_fp:
		config = ConfigObj(args.config)
		args.output_fp = config["SNP_DATABASE_FP"]
	
	build_snp_index(args.snp_dir,args.id_col,args.chr_col,args.locus_col,args.output_fp)
	
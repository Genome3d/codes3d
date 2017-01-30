from ftplib import FTP
import argparse, gzip, os, re, requests, subprocess, sys, wget

def download_snp_data():
	local_path = os.path.dirname(os.path.realpath(__file__)) + "/../lib/snps"

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	to_dl = "bed_chr_([0-9]+|[XY]|MT).*"
	ftp = FTP('ftp.ncbi.nih.gov')
	ftp.login()
	ftp.cwd('snp/organisms/human_9606_b146_GRCh37p13/BED')

	files = ftp.nlst()
	for file in files:
		if re.match(to_dl,file) and re.match(to_dl,file).group(0):
			print "Retrieving " + file
			with open("%s/%s" % (local_path, file), 'w') as local_file:
				ftp.retrbinary('RETR ' + file,local_file.write)

	ftp.quit()

def download_hic_data(to_dl):
	cell_lines = {}
	cell_lines["GM12878"] = ["GSM1551552"]
	cell_lines["HMEC"] = []
	cell_lines["HUVEC"] = []
	cell_lines["IMR90"] = []
	cell_lines["K562"] = []
	cell_lines["KBM7"] = []
	cell_lines["NHEK"] = []

	#hic_dir = os.path.dirname(os.path.realpath(__file__)) + "/../lib/hic_data"
	hic_dir = "/media/cam/UUI/lib/hic_data"

	ftp = FTP('ftp.ncbi.nih.gov')
	ftp.login()
	ftp.cwd('geo/samples/GSM1551nnn')
	for key in cell_lines.keys():
		if (to_dl and key in to_dl) or not to_dl:
			local_path = '%s/%s' % (hic_dir,key)
			if not os.path.isdir(local_path):
				os.makedirs(local_path)
			for dataset in cell_lines[key]:
				ftp.cwd(dataset + '/suppl')
				regex = dataset + "_HIC\d{3}_merged_nodups.txt.gz"
				files = ftp.nlst()
				for file in files:
					if re.match(regex,file):
						print "Retrieving " + file
						with open("%s/%s" % (local_path, file), 'w') as local_file:
							ftp.retrbinary('RETR ' + file,local_file.write)
	ftp.quit()

def download_gene_reference():
	local_path = os.path.dirname(os.path.realpath(__file__)) + "/../lib"

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	res = requests.get('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/reference/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz')
	with open(local_path + '/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz','wb') as write_gzip:
		write_gzip.write(res.content)

	with gzip.open(local_path + '/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz','rb') as gzip_in:
		genes_file = gzip_in.read()

	with open(local_path + '/gencode.v19.genes.v6p_model.patched_contigs.gtf','w') as genes_out:
		genes_out.write(genes_file)

	os.remove(local_path + '/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz')

def download_cis_eqtls():
	local_path = os.path.dirname(os.path.realpath(__file__)) + "/../lib/eQTLs"

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	eqtl_url = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/single_tissue_eqtl_data/GTEx_Analysis_V6_eQTLs.tar.gz"

	#subprocess.call(["wget",eqtl_url,"-P",local_path])
	subprocess.call(["tar","-xzf",local_path + "/GTEx_Analysis_V6_eQTLs.tar.gz","-C",local_path])

def download_human_genome():
	#local_path = os.path.dirname(os.path.realpath(__file__)) + "/../lib/hg19"
	local_path = "/media/cam/UUI/lib/hg19"

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	to_dl = "Homo_sapiens\.GRCh37\.75\.dna\.chromosome\.([0-9]+|[XY]|MT)\.fa\.gz"
	ftp = FTP('ftp.ensembl.org')
	ftp.login()
	ftp.cwd('pub/release-75/fasta/homo_sapiens/dna')

	files = ftp.nlst()
	for file in files:
		if re.match(to_dl,file) and re.match(to_dl,file).group(0):
			print "Retrieving " + file
			with open("%s/%s" % (local_path, file), 'w') as local_file:
				ftp.retrbinary('RETR ' + file,local_file.write)

	ftp.quit()

def download_expression_data():
	local_path = "/media/cam/UUI/lib"

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	eqtl_url = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"

	subprocess.call(["wget",eqtl_url,"-P",local_path])


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Download the default data files for CoDeS3D")
	parser.add_argument("-s","--snps",action="store_true",default=False,help="Download dbSNP data (build 146).")
	parser.add_argument("-c","--hic",nargs='*',help="Download Hi-C data. If no arguments are given, this will download Hi-C datasets for cell-lines\
		GM12878, HMEC, HUVEC, IMR90, K562, KBM7, and NHEK. Additionally, any of these can be passed as an argument in a space-separated list, e.g.\
		`-c GM12878 NHEK`. NOTE: THIS COMMAND DOWNLOADS TERABYTES WORTH OF DATA, AND MAY RUN FOR A LONG TIME. Consider downloading individual cell lines.")
	parser.add_argument("-g","--gene",action="store_true",default=False,help="Download GENCODE gene reference from GTEx portal.")
	parser.add_argument("-e","--cis_eqtls",action="store_true",default=False,help="Download cis-eQTLs found in GTEx analysis v6.")
	parser.add_argument("-m","--hg19",action="store_true",default=False,help="Download GRCh37.p13 (hg19) build of the human genome.")
	parser.add_argument("-x","--expression_data",action="store_true",default=False)
	args = parser.parse_args()

	if args.snps:
		download_snp_data()
	elif args.hic:
		download_hic_data(args.hic)
	elif args.gene:
		download_gene_reference()
	elif args.cis_eqtls:
		download_cis_eqtls()
	elif args.hg19:
		download_human_genome()
	elif args.expression_data:
		download_expression_data()
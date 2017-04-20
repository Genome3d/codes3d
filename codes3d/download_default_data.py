#!/usr/bin/env python

from ftplib import FTP
import argparse, codes3d, configparser, gzip, os, re, requests, subprocess, sys

def download_snp_data(conf, do_not_build_dbs):
	local_path = os.path.join(conf.get("Defaults","LIB_DIR"),"snps")

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
			with open(os.path.join(local_path, file), 'w') as local_file:
				ftp.retrbinary('RETR ' + file,local_file.write)
			print "Extracting " + file
			subprocess.call(["gzip","-d",os.path.join(local_path,file)])

	ftp.quit()
	if not do_not_build_dbs:
		codes3d.build_snp_index(local_path,os.path.join(conf.get("Defaults","LIB_DIR"),"snp_index_b146.db"),conf)

def download_hic_data(to_dl,conf,do_not_build_dbs):
	cell_lines = {}
	cell_lines["GM12878"] = ["GSM1551552","GSM1551553","GSM1551554","GSM1551555","GSM1551556","GSM1551557","GSM1551558","GSM1551559","GSM1551560","GSM1551561","GSM1551562","GSM1551563","GSM1551564","GSM1551565","GSM1551566","GSM1551567","GSM1551568","GSM1551569","GSM1551570","GSM1551571","GSM1551572","GSM1551573","GSM1551574"]
	cell_lines["HeLa"] = ["GSM1551632","GSM1551633","GSM1551634","GSM1551635"]
	cell_lines["HMEC"] = ["GSM1551607","GSM1551608","GSM1551609","GSM1551610","GSM1551611","GSM1551612"]
	cell_lines["HUVEC"] = ["GSM1551629","GSM1551630","GSM1551631"]
	cell_lines["IMR90"] = ["GSM1551599","GSM1551600","GSM1551601","GSM1551602","GSM1551603","GSM1551604","GSM1551605"]
	cell_lines["K562"] = ["GSM1551618","GSM1551619","GSM1551620","GSM1551621","GSM1551622","GSM1551623"]
	cell_lines["KBM7"] = ["GSM1551624","GSM1551625","GSM1551626","GSM1551627","GSM1551628"]
	cell_lines["NHEK"] = ["GSM1551614"]#,"GSM1551615","GSM1551616"]

	hic_dir = conf.get("Defaults","HIC_DATA_DIR")

	filelist = {}

	print "Retrieving file list..."
	ftp_domain = 'ftp.ncbi.nih.gov'
	ftp_dir = 'geo/samples/GSM1551nnn'
	ftp = FTP(ftp_domain)
	ftp.login()
	ftp.cwd(ftp_dir)
	for key in cell_lines.keys():
		if (to_dl and key in to_dl) or not to_dl:
			filelist[key] = []
			for dataset in cell_lines[key]:
				ftp.cwd(dataset + '/suppl')
				regex = dataset + "_HIC\d{3}_merged_nodups.txt.gz"
				files = ftp.nlst()
				for file in files:
					if re.match(regex,file):
						filelist[key].append((dataset,file))
				ftp.cwd('../../')
	ftp.quit()

	print "Downloading files..."
	for cell_line, datasets in filelist.items():
		local_path = os.path.join(hic_dir,cell_line)
		if not os.path.isdir(local_path):
			os.makedirs(local_path)
		for dataset in datasets:
			print "\tRetrieving " + dataset[0]
			fileurl = "ftp://%s/%s/%s/suppl/%s" % (ftp_domain,ftp_dir,dataset[0],dataset[1])
			subprocess.call(["wget",fileurl,"-P",local_path])

	print "Processing files..."
	for cell_line, datasets in filelist.items():
		local_path = os.path.join(hic_dir,cell_line)
		for dataset in datasets:
			print "\tExtracting " + dataset[1]
			subprocess.call(["gzip","-dk",os.path.join(local_path,dataset[1])])
			if not do_not_build_dbs:
				codes3d.build_hic_index(os.path.join(local_path,dataset[1][:dataset[1].rfind('.gz')]))


def download_gene_reference(conf,do_not_build_dbs):
	local_path = local_path = conf.get("Defaults","LIB_DIR")

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	res = requests.get('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/reference/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz')
	gene_fp = os.path.join(local_path,'gencode.v19.genes.v6p_model.patched_contigs.gtf')
	with open(gene_fp + ".gz",'wb') as write_gzip:
		write_gzip.write(res.content)

	with gzip.open(gene_fp + ".gz",'rb') as gzip_in:
		genes_file = gzip_in.read()

	with open(gene_fp,'w') as genes_out:
		genes_out.write(genes_file)

	os.remove(gene_fp + ".gz")
	if not do_not_build_dbs:
		codes3d.build_gene_index([gene_fp],gene_fp + ".bed",gene_fp + ".db",conf)

def download_cis_eqtls(conf,do_not_build_dbs):
	local_path = local_path = os.path.join(conf.get("Defaults","LIB_DIR"))

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	eqtl_url = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/single_tissue_eqtl_data/GTEx_Analysis_V6_eQTLs.tar.gz"

	subprocess.call(["wget",eqtl_url,"-P",local_path])
	subprocess.call(["tar","-xzf",os.path.join(local_path,"/GTEx_Analysis_V6_eQTLs.tar.gz"),"-C",local_path])
	if not do_not_build_dbs:
		for snpgene in os.listdir(os.path.join(local_path,"/eqtls")):
			codes3d.build_eqtl_index(snpgene)

def download_human_genome(conf,do_not_build_dbs):
	local_path = os.path.join(conf.get("Defaults","LIB_DIR"),"hg19")

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
			with open(os.path.join(local_path, file), 'w') as local_file:
				ftp.retrbinary('RETR ' + file,local_file.write)

	ftp.quit()

	print "Decompressing and concatenating files (please be patient, this may take some time)..."
	subprocess.call("zcat %s/*.gz > Homo_sapiens.GRCh37.75.dna.fa" % local_path,shell=True)

	if not do_not_build_dbs:
		codes3d.digest_genome(os.path.join(local_path,"Homo_sapiens.GRCh37.75.dna.fa"),"MboI",os.path.join(conf.get("Defaults","LIB_DIR"),"Homo_sapiens.GRCh37.75.dna.fragments.bed"),os.path.join(conf.get("Defaults","LIB_DIR"),"Homo_sapiens.GRCh37.75.dna.fragments.db"))


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Download the default data files for CoDeS3D. These will be downloaded into the default directories specified in the supplied config file (docs/conf.py by default).")
	parser.add_argument("-s","--snps",action="store_true",default=False,help="Download dbSNP data (build 146).")
	parser.add_argument("-i","--hic",nargs='*',help="Download Hi-C data. If no arguments are given, this will download Hi-C datasets for cell-lines\
		GM12878, HeLa, HMEC, HUVEC, IMR90, K562, KBM7, and NHEK. Additionally, any of these can be passed as an argument in a space-separated list, e.g.\
		`-i GM12878 NHEK`. NOTE: THIS COMMAND DOWNLOADS TERABYTES WORTH OF DATA, AND MAY RUN FOR A LONG TIME. Consider downloading individual cell lines.")
	parser.add_argument("-g","--gene",action="store_true",default=False,help="Download GENCODE gene reference from GTEx portal.")
	parser.add_argument("-e","--cis_eqtls",action="store_true",default=False,help="Download cis-eQTLs found in GTEx analysis v6.")
	parser.add_argument("-m","--hg19",action="store_true",default=False,help="Download GRCh37.p13 (hg19) build of the human genome.")
	parser.add_argument("-x","--expression_data",action="store_true",default=False)
	parser.add_argument("-c","--config_file",default=os.path.join(os.path.dirname(__file__),"../docs/codes3d.conf"),help="The configuration file to use to resolve library directories.")
	parser.add_argument("-b","--do_not_build_dbs",action="store_true",default=False,help="Do not build associated databases/indices for downloaded data (default: False).")
	args = parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config_file)

	if args.snps:
		download_snp_data(config,args.do_not_build_dbs)
	if args.hic:
		download_hic_data(args.hic,config,args.do_not_build_dbs)
	if args.gene:
		download_gene_reference(config,args.do_not_build_dbs)
	if args.cis_eqtls:
		download_cis_eqtls(config,args.do_not_build_dbs)
	if args.hg19:
		download_human_genome(config,args.do_not_build_dbs)

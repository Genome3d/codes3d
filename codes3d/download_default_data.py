#!/usr/bin/env python

from ftplib import FTP
import argparse, codes3d, configparser, gzip, os, re, requests, subprocess, sys

# Tested, working
def download_snp_data(conf, do_not_build_dbs):
	local_path = os.path.join(conf.get("Defaults","LIB_DIR"),"snps")

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	nih_domain = "ftp.ncbi.nih.gov"
	ftp_dir = "snp/organisms/human_9606_b146_GRCh37p13/BED"
	to_dl = "bed_chr_([0-9]+|[XY]|MT).*"
	file_list = []
	ftp = FTP(nih_domain)
	ftp.login()
	ftp.cwd(ftp_dir)
	files = ftp.nlst()
	for file in files:
		if re.match(to_dl,file) and re.match(to_dl,file).group(0):
			file_list.append(file)
	ftp.quit()
	for file in file_list:
		print "Retrieving " + file
		snp_url = "ftp://%s/%s/%s" % (nih_domain,ftp_dir,file)
		subprocess.call(["wget",snp_url,"-P",local_path])
		print "Extracting " + file
		subprocess.call(["gzip","-d",os.path.join(local_path,file)])
	if not do_not_build_dbs:
		codes3d.build_snp_index(local_path,os.path.join(conf.get("Defaults","LIB_DIR"),"snp_index_dbSNP_b146.db"),conf)

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
	nih_domain = "ftp.ncbi.nih.gov"
	ftp_dir = "geo/samples/GSM1551nnn"

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
		for dataset_file in datasets:
			print "\tRetrieving " + dataset_file[0]
			fileurl = "ftp://%s/%s/%s/suppl/%s" % (ftp_domain,ftp_dir,dataset_file[0],dataset_file[1])
			subprocess.call(["wget",fileurl,"-P",local_path])

	print "Processing files..."
	for cell_line, datasets in filelist.items():
		local_path = os.path.join(hic_dir,cell_line)
		for dataset_file in datasets:
			print "\tExtracting " + dataset_file[1]
			subprocess.call(["gzip","-dk",os.path.join(local_path,dataset_file[1])])
			if not do_not_build_dbs:
				codes3d.build_hic_index(os.path.join(local_path,dataset_file[1][:dataset_file[1].rfind('.gz')]))

#Tested, working
def download_gene_reference(conf,do_not_build_dbs):
	local_path = local_path = conf.get("Defaults","LIB_DIR")

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	gene_ref_url = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/reference/gencode.v19.genes.v6p_model.patched_contigs.gtf.gz"
	print "Retrieving GTEx gene reference..."
	subprocess.call(["wget",gene_ref_url,"-P",local_path])
	print "Extracting file..."
	gene_fp = os.path.join(local_path,'gencode.v19.genes.v6p_model.patched_contigs.gtf')
	subprocess.call(["gzip","-dk",gene_fp+".gz"])
	if not do_not_build_dbs:
		codes3d.build_gene_index([gene_fp],gene_fp + ".bed",gene_fp + ".db",conf)
		print "Cleaning up..."
		os.remove(gene_fp + ".gz")
		os.remove(gene_fp)

#Tested, working
def download_cis_eqtls(conf,do_not_build_dbs):
	lib_dir = conf.get("Defaults","LIB_DIR")
	eqtl_dir = conf.get("Defaults","EQTL_DATA_DIR")
	if not os.path.isdir(eqtl_dir):
		os.makedirs(eqtl_dir)

	eqtl_url = "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/single_tissue_eqtl_data/GTEx_Analysis_V6_eQTLs.tar.gz"

	print "Retrieving eQTL file..."
	subprocess.call(["wget",eqtl_url,"-P",lib_dir])
	print "Extracting GTEx_Analysis_V6_eQTLs.tar.gz..."
	subprocess.call(["tar","-xzf",os.path.join(lib_dir,"GTEx_Analysis_V6_eQTLs.tar.gz"),"-C",eqtl_dir])
	if not do_not_build_dbs:
		print "Files extracted, proceeding to build databases."
		for snpgene in os.listdir(eqtl_dir):
			if snpgene.endswith(".snpgene"):
				print "\tBuilding %s.db..." % snpgene[:snpgene.rfind('.')];
				codes3d.build_eqtl_index(os.path.join(eqtl_dir,snpgene));
				print "Cleaning up..."
				for snpgene in os.listdir(eqtl_dir):
					if snpgene.endswith(".snpgene"):
						os.remove(os.path.join(eqtl_dir,snpgene))
				os.remove(os.path.join(lib_dir,"GTEx_Analysis_V6_eQTLs.tar.gz"))

# Tested, working
def download_human_genome(conf,do_not_build_dbs):
	local_path = os.path.join(conf.get("Defaults","LIB_DIR"),"hg19")

	if not os.path.isdir(local_path):
		os.makedirs(local_path)

	ensembl_domain = "ftp.ensembl.org"
	ftp_dir = "pub/release-75/fasta/homo_sapiens/dna"
	file_list = []
	to_dl = "Homo_sapiens\.GRCh37\.75\.dna\.chromosome\.([0-9]+|[XY]|MT)\.fa\.gz"
	ftp = FTP(ensembl_domain)
	ftp.login()
	ftp.cwd(ftp_dir)

	files = ftp.nlst()
	for file in files:
		if re.match(to_dl,file) and re.match(to_dl,file).group(0):
			file_list.append(file)
	ftp.quit()
	for file in file_list:
		print "Retrieving " + file
		chr_url = "ftp://%s/%s/%s" % (ensembl_domain,ftp_dir,file)
		subprocess.call(["wget",chr_url,"-P",local_path])

	print "Decompressing and concatenating files (please be patient, this may take some time)..."
	subprocess.call("zcat %s/*.gz > %s/Homo_sapiens.GRCh37.75.dna.fa" % local_path,shell=True)

	if not do_not_build_dbs:
		codes3d.digest_genome(os.path.join(local_path,"Homo_sapiens.GRCh37.75.dna.fa"),"MboI",os.path.join(conf.get("Defaults","LIB_DIR"),"Homo_sapiens.GRCh37.75.dna.fragments.bed"),os.path.join(conf.get("Defaults","LIB_DIR"),"Homo_sapiens.GRCh37.75.dna.fragments.db"))


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Download the default data files for CoDeS3D. These will be downloaded into the default directories specified in the supplied config file (docs/conf.py by default).")
	parser.add_argument("-s","--snps",action="store_true",default=False,help="Download dbSNP data (build 146).")
	parser.add_argument("-i","--hic",nargs='*',help="Download Hi-C data. If no arguments are given, this will download Hi-C datasets for cell-lines\
		GM12878, HeLa, HMEC, HUVEC, IMR90, K562, KBM7, and NHEK. Additionally, any of these can be passed as an argument in a space-separated list, e.g.\
		`-i GM12878 NHEK`. NOTE: THIS COMMAND HUNDREDS OF GIGABYTES OF DATA, AND MAY RUN FOR A LONG TIME. Consider downloading individual cell lines.")
	parser.add_argument("-g","--gene",action="store_true",default=False,help="Download GENCODE gene reference from GTEx portal.")
	parser.add_argument("-e","--cis_eqtls",action="store_true",default=False,help="Download cis-eQTLs found in GTEx analysis v6.")
	parser.add_argument("-p","--hg19",action="store_true",default=False,help="Download GRCh37.p13 (hg19) build of the human genome.")
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

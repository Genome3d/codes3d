#!/usr/bin/env python

import argparse
import ast
import codes3d
import configparser
import os
import multiprocessing
import sqlite3
import time

def parse_failed_requests(requests_fp):
	reqLists = []
	with open(requests_fp,'r') as requests:
		for line in requests:
			reqLists.append(ast.literal_eval(line))
	return reqLists

def get_GTEx_responses(reqLists,snps,genes,gene_database_fp,num_processes,input_dir):
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	manager = multiprocessing.Manager()
	eqtls = {}
	print("\t\tRequests to send: " + str(len(reqLists)))
	gtexResponses = manager.list()
	procPool = multiprocessing.Pool(processes=num_processes)
	for i,reqList in enumerate(reqLists,start=1):
		print("\t\t\tSending request %s of %s" % (i,len(reqLists)))
		procPool.apply_async(codes3d.send_GTEx_query, (reqList,gtexResponses))
		time.sleep(10)
	procPool.close()
	procPool.join()
	print("\t\tNumber of GTEx responses received: " + str(len(gtexResponses)))
	results = []
	failed_requests = []
	for response in gtexResponses:
		try:
			results += response[1].json()["result"]
		except Exception as e:
			print("\t\tWARNING: bad response.\n\t\t\tException: %s\n\t\t\tResponse:\n\t\t\t%s" % (e,response[1]))
			failed_requests.append(response[0])
	if failed_requests:
		with open(input_dir + "/failed_GTEx_requests.txt",'w') as failed_requests_file:
			failed_requests_file.write(str(failed_requests))
	for result in results:
		geneSymbol = result["geneSymbol"]
		snpId = result["variantId"]
		if str(geneSymbol) == "gene not found" or not snps.has_key(snpId):
			continue
		if not eqtls.has_key(snpId):
			eqtls[snpId] = {}
			eqtls[snpId]["snp_info"] = { "chr": snps[snpId]["chr"], "locus": snps[snpId]["locus"] }
			num_sig[snpId] = 0
		if not eqtls[snpId].has_key(geneSymbol):
			gene_chr = "NA"
			gene_start = "NA"
			gene_end = "NA"
			p_thresh = "NA"
			
			max_length = 0
			for gene_stat in gene_index.execute("SELECT chr, start, end, p_thresh FROM genes WHERE symbol=?", [geneSymbol]):
				if gene_stat[2] - gene_stat[1] > max_length: #Consider "canonical" to be the longest record where multiple records are present
					if gene_stat[0].startswith("chr"):
						gene_chr = gene_stat[0][gene_stat[0].find("chr")+3:]
					else:
						gene_chr = gene_stat[0]
					gene_start = gene_stat[1]
					gene_end = gene_stat[2]
					p_thresh = gene_stat[3]
					max_length = gene_stat[2] - gene_stat[1]
			cis = gene_chr == eqtls[snpId]["snp_info"]["chr"] and (eqtls[snpId]["snp_info"]["locus"] > gene_start - 1000000 and eqtls[snpId]["snp_info"]["locus"] < gene_end + 1000000) #eQTL is cis if the SNP is within 1Mbp of the gene
			
			eqtls[snpId][geneSymbol] = {}
			eqtls[snpId][geneSymbol]["ens_id"] = result["gencodeId"]
			eqtls[snpId][geneSymbol]["gene_chr"] = gene_chr
			eqtls[snpId][geneSymbol]["gene_start"] = gene_start
			eqtls[snpId][geneSymbol]["gene_end"] = gene_end
			eqtls[snpId][geneSymbol]["cell_lines"] = list(genes[snpId][geneSymbol])
			eqtls[snpId][geneSymbol]["p_thresh"] = p_thresh
			eqtls[snpId][geneSymbol]["tissues"] = {}
			if cis:
				eqtls[snpId][geneSymbol]["cis?"] = True
			else:
				eqtls[snpId][geneSymbol]["cis?"] = False
		if not result["pvalue"] == "NA":
			p = float(result["pvalue"])
			eqtls[snpId][geneSymbol]["tissues"][result["tissueId"]] = {"pvalue": p}
	
	print("Writing to eQTL file...")
	if os.path.isfile(input_dir + "/eqtls.txt"):
		upsert = input("Appending output to existing eQTL file (%s). Continue? [y/N] " % input_dir + "/eqtls.txt")
		if not upsert.lower() == 'y':
			print("Did not write to existing SNP database.")
			return
	with open(input_dir + "/eqtls.txt",'a') as eqtlfile:
		for snp in eqtls.keys():
			for gene in eqtls[snp].keys():
				if not gene == "snp_info":
					for tissue in eqtls[snp][gene]["tissues"].keys():
						eqtlfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snp,eqtls[snp]["snp_info"]["chr"],eqtls[snp]["snp_info"]["locus"],gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],eqtls[snp][gene]["cell_lines"],eqtls[snp][gene]["cis?"],eqtls[snp][gene]["p_thresh"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"]))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Tool for re-sending failed GTEx requests, and optionally appending the results to an existing eQTLs file.")
	parser.add_argument("-i","--input_dir",required=True,help="The output directory from the CoDeS3D run containing failed_requests.")
	parser.add_argument("-c","--config",default="docs/codes3d.conf.py",help="The configuration file to be used in this instance (default: conf.py)")
	parser.add_argument("-p","--num_processes",type=int,default=1,help="Desired number of processes for multiprocessing (default: 1).")
	args = parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config)
	gene_database_fp = config.get("Defaults","GENE_DATABASE_FP")
	reqLists = parse_failed_requests(args.input_dir + "/failed_GTEx_requests.txt")
	snps = codes3d.parse_snps_files(args.input_dir + "/snps.txt")
	genes = codes3d.parse_genes_files(args.input_dir + "/genes.txt")
	get_GTEx_responses(reqLists,snps,genes,gene_database_fp,args.num_processes,args.input_dir)

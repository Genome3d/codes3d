#!/usr/bin/python
import ast
import csv
import argparse
import os
import sqlite3
import shelve
import pybedtools
import requests
import multiprocessing
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.ticker import FuncFormatter
from itertools import cycle
import re
import pandas
import json
import bisect

from configobj import ConfigObj
from sets import Set

from wikipathways_api_client import WikipathwaysApiClient

def process_inputs(inputs,snp_database_fp,output_dir,suppress_intermediate_files=False):
	print "Processing input..."
	snp_db = sqlite3.connect(snp_database_fp)
	snp_db.text_factory = str
	snp_index = snp_db.cursor()
	snps = {}
	for input in inputs:
		if os.path.isfile(input):
			with open(input,'r') as infile:
				for line in infile:
					id = line.strip().split(' ')[0]
					snp = None
					if id.startswith("rs"):
						snp_index.execute("SELECT * FROM snps WHERE rsID=?",(id,))
					else:
						chr = id[id.find("chr")+3:id.find(':')]
						locus = int(id[id.find(':')+1:])
						snp_index.execute("SELECT * FROM snps WHERE chr=? and locus=?",[chr,locus])
					snp = snp_index.fetchone()
					if snp == None:
						print "Warning: %s does not exist in SNP database." % id
					else:
						snps[snp[0]]=(snp[1],snp[2])
		else:
			snp = None
			if input.startswith("rs"):
				snp_index.execute("SELECT * FROM snps WHERE rsID=?",(input,))
			else:
				chr = input[input.find("chr")+3:input.find(':')]
				locus = int(input[input.find(':')+1:])
				snp_index.execute("SELECT * FROM snps WHERE chr=? and locus=?",[chr,locus])
			snp = snp_index.fetchone()
			if snp == None:
				print "Warning: %s does not exist in SNP database." % input
			else:
				snps[snp[0]]=(snp[1],snp[2])
	if not suppress_intermediate_files:
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		with open(output_dir + "/snps.txt",'w') as snpfile:
			for snp,loc in snps.items():
				snpfile.write("%s\t%s\t%s\n" % (snp,loc[0],loc[1]))
	return snps

def find_interactions(snps,fragment_database_fp,hic_data_dir,include,exclude,output_dir,suppress_intermediate_files=False):
	print "Finding interactions..."
	#Query fragmentIndex to find to which fragment the SNP belongs
	fragment_index_db = sqlite3.connect(fragment_database_fp)
	fragment_index_db.text_factory = str
	fragment_index = fragment_index_db.cursor()
	#Look for all interactions involving said fragment in the HiC databases
	if include:
		include_set = Set(include)
	if exclude:
		exclude_set = Set(exclude)
	interactions = {} #A mapping of each SNP to the fragments with which the fragment it is on interacts
	for snp in snps.keys():
		interactions[snp] = {}
	for cell_line in os.listdir(hic_data_dir):
		if (include and not cell_line in include) or (exclude and cell_line in exclude):
			continue
		if os.path.isdir(hic_data_dir + '/' + cell_line):
			interactions[snp][cell_line] = Set([])
			print "\tSearching cell line " + cell_line
			for replicate in os.listdir(hic_data_dir + '/' + cell_line):
				if replicate.endswith(".db"):
					rep_ints = shelve.open(hic_data_dir + '/' + cell_line + '/' + replicate)
					print "\t\tSearching replicate " + replicate
					for snp in snps.keys():
                                                interactions[snp][cell_line] = Set([])
						print "\t\t\tFinding interactions for " + str(snp)
						fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",["chr" + snps[snp][0],snps[snp][1],snps[snp][1]])
						snp_fragment_result = fragment_index.fetchone()
						if snp_fragment_result == None:
							print "Warning: error retrieving SNP fragment for SNP " + snp
							continue
						snp_fragment = snp_fragment_result[0]
						snp_chr = snps[snp][0]
						key = "%s.%s" % (snp_chr,snp_fragment)
						frags = rep_ints[key]
						for interaction in frags:
							interactions[snp][cell_line].add(interaction)
	if not suppress_intermediate_files:
		with open(output_dir + "/snp-gene_interactions.txt",'w') as intfile:
			for snp in interactions.keys():
				for cell_line in interactions[snp].keys():
					for interaction in interactions[snp][cell_line]:
						intfile.write("%s\t%s\t%s\t%s\n" % (snp,cell_line,interaction[0],interaction[1]))
	return interactions

def find_genes(interactions,fragment_database_fp,gene_bed_fp,output_dir,suppress_intermediate_files=False):
	print "Identifying interactions with genes..."
	fragment_index_db = sqlite3.connect(fragment_database_fp)
	fragment_index_db.text_factory = str
	fragment_index = fragment_index_db.cursor()
	hs_gene_bed = pybedtools.BedTool(gene_bed_fp)
	genes = {}
	for snp in interactions.keys():
		#Generate BED file of all fragments interacting with SNP-containing fragment
		for cell_line in interactions[snp].keys():
			snpgenes_exist = False
			temp_snp_bed = open(output_dir + "/temp_snp_bed.bed",'w')
			for interaction in interactions[snp][cell_line]:
				fragment_index.execute("SELECT start, end FROM fragments WHERE chr=? and fragment=?",["chr" + interaction[0],interaction[1]])
				fragment_pos = fragment_index.fetchone()
				if fragment_pos == None:
					print "\tWarning: error retrieving fragment %s on chromosome %s" % (interaction[1],interaction[0])
					continue
				temp_snp_bed.write("%s\t%s\t%s\n" % ("chr" + interaction[0],fragment_pos[0],fragment_pos[1]))
				if not snpgenes_exist:
					snpgenes_exist = True
			temp_snp_bed.close()
			if snpgenes_exist:
				if not genes.has_key(snp):
					genes[snp] = {}
				int_bed = pybedtools.BedTool(output_dir + "/temp_snp_bed.bed")
				#Get intersection of this BED file with BED file detailing gene locations
				gene_bed = hs_gene_bed.intersect(int_bed,u=True)
				#Return a list of genes with which SNP is interacting
				for feat in gene_bed:
					if not genes[snp].has_key(feat.name):
						genes[snp][feat.name] = Set([])
					genes[snp][feat.name].add(cell_line)
	os.remove(output_dir + "/temp_snp_bed.bed")
	snps_to_remove = []
	for snp in interactions.keys():
		if not genes.has_key(snp):
			print "\tNo SNP-gene spatial interactions detected for %s, removing from analysis" % (snp,)
			snps_to_remove.append(snp)
	for snp in snps_to_remove:
		del snps[snp]
		del interactions[snp]
	if not suppress_intermediate_files:
		with open(output_dir + "/genes.txt",'w') as genefile:
			for snp in genes.keys():
				for gene in genes[snp].keys():
					for cell_line in genes[snp][gene]:
						genefile.write("%s\t%s\t%s\n" % (snp,gene,cell_line))
	return genes

def find_eqtls(snps,genes,eqtl_data_dir,gene_database_fp,fdr_threshold,local_databases_only,num_processes,output_dir,suppress_intermediate_files=False):
	print "Identifying eQTLs of interest..."
	eqtls = {} #A mapping of SNPs to genes with which they have an eQTL relationship, which in turn maps to a list of tissues in which this eQTL occurs
	p_values = [] #A sorted list of all p-values for use computing FDR
	num_tests = 0 #Total number of tests done
	num_sig = {} #Number of eQTLs deemed significant under the given threshold
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	print "\tQuerying local databases."
	for db in os.listdir(eqtl_data_dir): #Iterate through databases of eQTLs by tissue type
		tissue = db[:db.rfind('.')]
		print "\t\tQuerying " + tissue
		eqtl_index_db = sqlite3.connect(eqtl_data_dir + '/' + db)
		eqtl_index_db.text_factory = str
		eqtl_index = eqtl_index_db.cursor()
		for snp in genes.keys():
			for gene in genes[snp].keys():
				num_tests += 1
				for eqtl in eqtl_index.execute("SELECT pvalue FROM eqtls WHERE rsID=? AND gene_symbol=?",(snp,gene)): #Pull down all eQTLs related to a given SNP to test for relevance:
					if not eqtls.has_key(snp):
						eqtls[snp] = {} #Stores eQTLs relevant to this SNP, and, for each eQTL, a list of the tissues in which it is found and the associated p-value with the eQTL in this tissue
						num_sig[snp] = 0
					if not eqtls[snp].has_key(gene):
						gene_chr = "NA"
						gene_start = "NA"
						gene_end = "NA"
						p_thresh = "NA"
						max_length = 0
						for gene_stat in gene_index.execute("SELECT chr, start, end, p_thresh FROM genes WHERE symbol=?", [gene]):
							if gene_stat[2] - gene_stat[1] > max_length: #Consider "canonical" to be the longest record where multiple records are present
								if gene_stat[0].startswith("chr"):
									gene_chr = gene_stat[0][gene_stat[0].find("chr")+3:]
								else:
									gene_chr = gene_stat[0]
								gene_start = gene_stat[1]
								gene_end = gene_stat[2]
								p_thresh = gene_stat[3]
								max_length = gene_stat[2] - gene_stat[1]
						if gene_chr == "NA":
							print "\t\tWarning: no entry in gene database for " + gene
						eqtls[snp][gene] = {}
						eqtls[snp][gene]["gene_chr"] = gene_chr
						eqtls[snp][gene]["gene_start"] = gene_start
						eqtls[snp][gene]["gene_end"] = gene_end
						eqtls[snp][gene]["cell_lines"] = list(genes[snp][gene])
						eqtls[snp][gene]["p_thresh"] = p_thresh
						eqtls[snp][gene]["tissues"] = {}
						cis = gene_chr == snps[snp][0] and (snps[snp][1] > gene_start - 1000000 and snps[snp][1] < gene_end + 1000000) #eQTL is cis if the SNP is within 1Mbp of the gene
						if cis:
							eqtls[snp][gene]["cis?"] = True
						else:
							eqtls[snp][gene]["cis?"] = False
					try:
						p = eqtl[0]
						eqtls[snp][gene]["tissues"][tissue] = {"pvalue": p}
						bisect.insort(p_values,p)
					except sqlite3.OperationalError:
						pass

	if not local_databases_only:
		print "\tQuerying GTEx online database."
		num_tests = get_GTEx_response(snps,genes,num_tests,num_sig,gene_database_fp,eqtls,p_values,num_processes,output_dir)

	snps_to_remove = []
	for snp in snps:
		if not eqtls.has_key(snp):
			print "\t\tWarning: no eQTLs found for SNP %s, removing from analysis." % snp
			del genes[snp]
			del interactions[snp]
			snps_to_remove.append(snp)
	for snp in snps_to_remove:
		del snps[snp]
	
	eqtlfile = None
	if not suppress_intermediate_files:
		print "\tComputing q-values and writing eqtlfile..."
		eqtlfile = open(output_dir + '/eqtls.txt','w')
	else:
		print "\tComputing q-values..."
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			for tissue in eqtls[snp][gene]["tissues"].keys():
				eqtls[snp][gene]["tissues"][tissue]["qvalue"] = compute_fdr(eqtls[snp][gene]["tissues"][tissue]["pvalue"],p_values,num_tests)
				if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
					num_sig[snp] += 1
				if not suppress_intermediate_files:
					eqtlfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snp,gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],eqtls[snp][gene]["cell_lines"],eqtls[snp][gene]["cis?"],eqtls[snp][gene]["p_thresh"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"],eqtls[snp][gene]["tissues"][tissue]["qvalue"]))
	
	if not suppress_intermediate_files:
		with open(output_dir + "/significant_eqtls_per_snp.txt",'w') as sig_eqtls:
			for snp,n in num_sig.items():
				sig_eqtls.write("%s\t%s\n" % (snp,n))
	
	return (eqtls,num_sig)

def get_GTEx_response(snps,genes,num_tests,num_sig,gene_database_fp,eqtls,p_values,num_processes,output_dir):
	tissues = Set(["Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Transformed_fibroblasts","Colon_Sigmoid","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood"])
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	manager = multiprocessing.Manager()
	reqLists = [[]]
	for snp in genes.keys():
		for gene in genes[snp].keys():
			for tissue in tissues:
				#We aren't interested in eQTLs already discovered from local databases
				if (not eqtls.has_key(snp) or (eqtls.has_key(snp) and not eqtls[snp].has_key(gene))) or (eqtls[snp].has_key(gene) and not eqtls[snp][gene]["tissues"].has_key(tissue)):
					if len(reqLists[-1]) < 1000:
						reqLists[-1].append({"snpId":snp,"gencodeId":gene,"tissueName":tissue})
					else:
						reqLists.append([{"snpId":snp,"gencodeId":gene,"tissueName":tissue}])
	print "\t\tRequests to send: " + str(len(reqLists))
	gtexResponses = manager.list()
	procPool = multiprocessing.Pool(processes=num_processes)
	for reqList in reqLists:
		procPool.apply_async(send_GTEx_query, (reqList,gtexResponses))
	procPool.close()
	procPool.join()
	print "\t\tNumber of GTEx responses received: " + str(len(gtexResponses))
	results = []
	failed_requests = open(output_dir + "/failed_requests.txt",'w')
	for response in gtexResponses:
		try:
			results += response[1].json()["result"]
		except Exception as e:
			print "\t\tWarning: bad response. Exception: %s" % e
			failed_requests.write(response[0])
	for result in results:
		geneSymbol = result["geneSymbol"]
		snpId = result["snpId"]
		if str(geneSymbol) == "gene not found" or not snps.has_key(snpId):
			continue
		num_tests += 1
		if not eqtls.has_key(snpId):
			eqtls[snpId] = {}
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
			cis = gene_chr == snps[snp][0] and (snps[snp][1] > gene_start - 1000000 and snps[snp][1] < gene_end + 1000000) #eQTL is cis if the SNP is within 1Mbp of the gene
			
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
			bisect.insort(p_values,p)
		
	return num_tests

def send_GTEx_query(reqList,gtexResponses):
	print "\t\tSending request to GTEx API..."
	try:
		gtexResponses.append((reqList,requests.post("http://gtexportal.org/api/v6/dyneqtl?v=clversion", json=reqList)))
	except requests.exceptions.ConnectionError:
		try:
			print "\t\tWarning: a connection error occurred with the GTEx service. Retrying..."
			gtexResponses.append((reqList,requests.post("http://gtexportal.org/api/v6/dyneqtl?v=clversion", json=reqList))) #Allow to crash if fails a second time
		except requests.exceptions.ConnectionError:
			print "\t\tRetry failed. Continuing, but results will be incomplete."
			gtexResponses.append((reqList,"Connection failure"))
			return
	print "\t\tResponse received."

def get_gene_expression_info(eqtls,expression_table_fp,output_dir,suppress_intermediate_files=False):
	print "Getting gene expression information..."
	gene_df = pandas.read_table(expression_table_fp,index_col='Symbol')
	gene_exp = {}
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			gene_exp[gene] = {}
			try:
				for tissue in list(gene_df.columns.values):
					gene_exp[gene][tissue] = gene_df.at[gene,tissue]
					gene_exp[gene]["max"] = gene_df.ix[gene].idxmax()
					if not isinstance(gene_exp[gene]["max"], str): #If the column header for max expression is not a string, the gene has multiple entries in the gene expression table
						gene_exp[gene]["max"] = gene_df.ix[gene].max().idxmax()
					gene_exp[gene]["min"] = gene_df.ix[gene].idxmin()
					if not isinstance(gene_exp[gene]["min"], str):
						gene_exp[gene]["min"] = gene_df.ix[gene].min().idxmin()
			except KeyError:
				print "Warning: Gene expression data not found for " + gene
				for tissue in list(gene_df.columns.values):
					gene_exp[gene][tissue] = "NA"
					gene_exp[gene]["max"] = "NA"
					gene_exp[gene]["min"] = "NA"
	if not suppress_intermediate_files:
		with open(output_dir + '/gene_expression.txt','w') as gene_exp_file:
			for gene in gene_exp.keys():
				for tissue in gene_exp[gene].keys():
					gene_exp_file.write("%s\t%s\t%s\n" % (gene,tissue,gene_exp[gene][tissue]))
	return gene_exp

def compute_fdr(p,p_values,num_tests):
	positives = bisect.bisect_right(p_values, p)
	return float(num_tests) * p / float(positives)

def produce_output(snps,genes,eqtls,num_sig,gene_exp,output_dir,suppress_graphs):
	#TODO: Tidy up this method using pandas
	print "Producing output..."
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)
	stat_table = open(output_dir + "/stat_table.txt",'w')
	stat_table.write("SNP\tChromosome\tLocus\tTotal_SNP-gene_Pairs\tTotal_eQTLs\n")
	summary = open(output_dir + "/summary.txt",'w')
	summary.write("SNP\tSNP_Chromosome\tSNP_Locus\tGene_Name\tGene_Chromosome\tGene_Start\tGene_End\tTissue\tp-value\tq-value\tCell_Lines\tGTEx_cis_p_Threshold\tcis_SNP-gene_interaction\tSNP-gene_Distance\tExpression_Level_In_eQTL_Tissue\tMax_Expressed_Tissue\tMaximum_Expression_Level\tMin_Expressed_Tissue\tMin_Expression_Level\n")
	for snp in snps.keys():
		stat_table.write(snp + '\t' + snps[snp][0] + '\t' + str(snps[snp][1]) + '\t' + str(len(genes[snp])) + '\t' + str(len(eqtls[snp])) + '\n')
		for gene in eqtls[snp].keys():
			distance_from_snp = 0
			if(not snps[snp][0] == eqtls[snp][gene]["gene_chr"]):
				distance_from_snp = "NA" #Not applicable to trans interactions
			elif(snps[snp][1] < eqtls[snp][gene]["gene_start"]):
				distance_from_snp = eqtls[snp][gene]["gene_start"] - snps[snp][1]
			elif(snps[snp][1] > eqtls[snp][gene]["gene_end"]):
				distance_from_snp = snps[snp][1] - eqtls[snp][gene]["gene_end"]
			
			gene_exp_max_tis = gene_exp[gene]["max"]
			if gene_exp_max_tis == "NA":
				gene_exp_max_val = "NA"
			else:
				gene_exp_max_val = gene_exp[gene][gene_exp_max_tis]
			gene_exp_min_tis = gene_exp[gene]["min"]
			if gene_exp_min_tis == "NA":
				gene_exp_min_val = "NA"
			else:
				gene_exp_min_val = gene_exp[gene][gene_exp_min_tis]
			eqtl_tissue = []
			for tissue in eqtls[snp][gene]["tissues"].keys():
				gene_exp_this_tissue = "NA"
				try:
					gene_exp_this_tissue = gene_exp[gene][tissue]
				except KeyError:
					print "\t\tWarning: No expression information for %s in %s" % (gene,tissue)
				eqtl_tissue.append((eqtls[snp][gene]["tissues"][tissue]["pvalue"],str(eqtls[snp][gene]["gene_chr"]),str(eqtls[snp][gene]["gene_start"]),str(eqtls[snp][gene]["gene_end"]),tissue,str(eqtls[snp][gene]["tissues"][tissue]["qvalue"]),eqtls[snp][gene]["cell_lines"],str(eqtls[snp][gene]["p_thresh"]),str(eqtls[snp][gene]["cis?"]),str(distance_from_snp),str(gene_exp_this_tissue)))
			eqtl_tissue.sort() #Sort by p-value
			for entry in eqtl_tissue:
				summary.write(snp + '\t' + snps[snp][0] + '\t' + str(snps[snp][1]) + '\t' + gene + '\t' + entry[1] + '\t' + entry[2] + '\t' + entry[3] + '\t' + entry[4] + '\t' + str(entry[0]) + '\t' + entry[5] + '\t')
				for i,cell_line in enumerate(entry[6],start=1):
					if not i == len(entry[6]):
						summary.write(cell_line + ',')
					else:
						summary.write(cell_line)
				summary.write('\t' + entry[7] +'\t' + entry[8] + '\t' + entry[9] + '\t' + entry[10] + '\t' + gene_exp_max_tis + '\t' + str(gene_exp_max_val) + '\t' + gene_exp_min_tis + '\t' + str(gene_exp_min_val) + '\n')
	stat_table.close()
	summary.close()
	if not suppress_graphs:
		produce_plots(snps,genes,eqtls,num_sig,output_dir)

def produce_plots(snps,genes,eqtls,num_sig,output_dir):
	print "\tProducing graphs"
	style.use("ggplot")
	if not os.path.isdir(output_dir + "/plots"):
		os.mkdir(output_dir + "/plots")
	int_colours = "rgb"
	eqtl_colours = "myc"
	int_colours = cycle(int_colours)
	eqtl_colours = cycle(eqtl_colours)
	snps_by_chr = {}
	for snp, info in snps.items():
		chrm = info[0]
		if not snps_by_chr.has_key(chrm):
			snps_by_chr[chrm] = []
		snps_by_chr[chrm].append((info[1],snp))
	
	int_colour_list = []
	eqtl_colour_list = []
	num_snpgenes = []
	num_eqtls = []
	rsIDs = []
	
	chrs = snps_by_chr.keys()
	chrs.sort(key=natural_keys) #So that the chromosomes are in a logical order on the graph
	chr_locs = []
	chr_ticks = []
	chrm_pos = 0
	last_count = 0
	count = 0
	
	for chrm in chrs:
		snp_list = snps_by_chr[chrm]
		snp_list.sort() #Sort by locus
		int_colour = int_colours.next()
		eqtl_colour = eqtl_colours.next()
		chr_locs.append(chrm_pos)
		chr_ticks.append(chrm)
		for snp in snp_list:
			num_snpgenes.append(len(genes[snp[1]]))
			num_eqtls.append(num_sig[snp[1]] * -1)
			rsIDs.append(snp[1])
			int_colour_list.append(int_colour)
			eqtl_colour_list.append(eqtl_colour)
			count += 1
		chrm_pos = chrm_pos + count
		count = 0
	
	plt.clf()
	plt.bar(range(len(num_snpgenes)),num_snpgenes,color=int_colour_list)
	plt.bar(range(len(num_eqtls)),num_eqtls,color=eqtl_colour_list)
	axes = plt.gca()
	axes.yaxis.set_major_formatter(FuncFormatter(abs_value_ticks)) #So that eQTL values won't appear negative on plot
	plt.vlines(chr_locs,axes.get_ylim()[0],axes.get_ylim()[1],colors="k")
	axes.set_xticks(range(len(rsIDs)))
	axes.set_xticklabels(rsIDs,rotation='vertical')
	ax2 = axes.twiny()
	ax2.set_xlim(axes.get_xlim())
	ax2.set_xticks(chr_locs)
	ax2.set_xticklabels(chr_ticks)
	plt.tight_layout()
	plt.savefig(output_dir + "/plots/snpgene_and_eqtl_overview.png",dpi=300,format="png")
	plt.clf()

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def abs_value_ticks(x, pos):
	return abs(x)

def interrogate_pathways(eqtls,fdr_threshold,output_dir):
	print "Retrieving pathway information..."
	pathways = {} #Keeps track of all pathways covered by genes with a statistically signficant eQTL relationship with a query SNP
	wikipathways = WikipathwaysApiClient()
	
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			for tissue in eqtls[snp][gene]["tissues"].keys():
				if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
					kwargs = {
						'query': gene,
						'organism': 'http://identifiers.org/taxonomy/9606' #Homo sapiens
					}
					pwresults = wikipathways.find_pathways_by_text(**kwargs)
					for pwres in pwresults:
						pwid = pwres["identifier"]
						if not pathways.has_key(pwid):
							pathways[pwid] = {}
							pathways[pwid]["name"] = pwres["name"]
							pathways[pwid]["genes"] = {}
						if not pathways[pwid]["genes"].has_key(gene):
							pathways[pwid]["genes"][gene] = {}
						if not pathways[pwid]["genes"][gene].has_key(snp):
							pathways[pwid]["genes"][gene][snp] = Set([])
						pathways[pwid]["genes"][gene][snp].add(tissue)

	with open(output_dir + "/pathways.txt",'w') as pwfile:
		pwfile.write("WikiPathways_ID\tPathway_Name\tGene_Symbol\tSNP\tTissue\n")
		for pwid in pathways.keys():
			for gene in pathways[pwid]["genes"].keys():
				for snp in pathways[pwid]["genes"][gene].keys():
					for tissue in pathways[pwid]["genes"][gene][snp]:
						pwfile.write("%s\t%s\t%s\t%s\t%s\n" % (pwid,pathways[pwid]["name"],gene,snp,tissue))
	return pathways

def parse_snps_file(snps_file):
	snps = {}
	with open(snps_file,'r') as snpfile:
		for line in snpfile:
			snp = line.strip().split('\t')
			snps[snp[0]] = (snp[1],int(snp[2]))
	return snps

def parse_interactions_file(interactions_file):
	interactions = {}
	with open(interactions_file,'r') as int_file:
		for line in int_file:
			interaction = line.strip().split('\t')
			snp = interaction[0]
			cell_line = interaction[1]
			if not interactions.has_key(snp):
				interactions[snp] = {}
			if not interactions[snp].has_key(cell_line):
				interactions[snp][cell_line] = Set([])
			interactions[snp][cell_line].add((interaction[2],int(interaction[3])))
	return interactions

def parse_genes_file(genes_file):
	genes = {}
	with open(genes_file,'r') as genefile:
		for line in genefile:
			gene = line.strip().split('\t')
			snp = gene[0]
			if not genes.has_key(snp):
				genes[snp] = {}
			if not genes[snp].has_key(gene[1]):
				genes[snp][gene[1]] = Set([])
			genes[snp][gene[1]].add(gene[2])
	return genes

def parse_eqtls_file(eqtl_file):
	eqtls = {}
	with open(eqtl_file,'r') as eqtlfile:
		for line in eqtlfile:
			eqtl = line.strip().split('\t')
			snp = eqtl[0]
			gene = eqtl[1]
			gene_chr = eqtl[2]
			try:
				gene_start = int(eqtl[3])
				gene_end = int(eqtl[4])
			except ValueError:
				gene_start = eqtl[3]
				gene_end = eqtl[4]
			cell_lines = ast.literal_eval(eqtl[5])
			cis = ast.literal_eval(eqtl[6])
			try:
				p_thresh = float(eqtl[7])
			except ValueError:
				p_thresh = eqtl[7]
			tissue = eqtl[8]
			p = float(eqtl[9])
			q = float(eqtl[10])
			if not eqtls.has_key(snp):
				eqtls[snp] = {}
			if not eqtls[snp].has_key(gene):
				eqtls[snp][gene] = {}
				eqtls[snp][gene]["gene_chr"] = gene_chr
				eqtls[snp][gene]["gene_start"] = gene_start
				eqtls[snp][gene]["gene_end"] = gene_end
				eqtls[snp][gene]["cell_lines"] = cell_lines
				eqtls[snp][gene]["p_thresh"] = p_thresh
				eqtls[snp][gene]["tissues"] = {}
				eqtls[snp][gene]["cis?"] = cis
			eqtls[snp][gene]["tissues"][tissue] = {"pvalue": p, "qvalue": q}
	return eqtls

def parse_number_significant_eqtls_per_snp_file(num_sig_file):
	num_sig = {}
	with open(num_sig_file,'r') as numsigfile:
		for line in numsigfile:
			ns = line.strip().split('\t')
			num_sig[ns[0]] = int(ns[1])
	return num_sig

def parse_gene_exp_file(gene_exp_file):
	gene_exp = {}
	with open(gene_exp_file,'r') as geneexpfile:
		for line in geneexpfile:
			ge = line.strip().split('\t')
			gene = ge[0]
			tissue = ge[1]
			try:
				exp = float(ge[2])
			except ValueError:
				exp = ge[2]
			if not gene_exp.has_key(gene):
				gene_exp[gene] = {}
			gene_exp[gene][tissue] = exp
	return gene_exp

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',required=True,help="The the dbSNP IDs or loci of SNPs of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-c","--config",default="conf.py",help="The configuration file to be used for this hiCquery run (default: conf.py)")
	parser.add_argument("-n","--include_cell_lines",nargs='+',help="Space-separated list of cell lines to include (others will be ignored). NOTE: Mutually exclusive with EXCLUDE_CELL_LINES.")
	parser.add_argument("-x","--exclude_cell_lines",nargs='+',help="Space-separated list of cell lines to exclude (others will be included). NOTE: Mutually exclusive with INCLUDE_CELL_LINES.")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	parser.add_argument("-l","--local_databases_only",action="store_true",default=False,help="Consider only local databases. Will only include cis-eQTLs if using downloadable GTEx dataset.")
	parser.add_argument("-s","--suppress_intermediate_files",action="store_true",default=False,help="Do not produce intermediate files. These can be used to run the pipeline from an intermediate stage in the event of interruption.")
	parser.add_argument("-p","--num_processes",type=int,default=1,help="Desired number of processes for multithreading (default: 1).")
	parser.add_argument("-g","--suppress_graphs",action="store_true",default=False,help="Do not produce plots of numbers of SNP-gene interactions versus numbers of significant eQTLs.")
	parser.add_argument("-f","--fdr_threshold",type=float,default=0.05,help="The FDR threshold to consider an eQTL statistically significant (default: 0.05).")
	args = parser.parse_args()
	config = ConfigObj(args.config)
	snp_database_fp = config["SNP_DATABASE_FP"]
	hic_data_dir = config["HIC_DATA_DIR"]
	fragment_bed_fp = config["FRAGMENT_BED_FP"]
	fragment_database_fp = config["FRAGMENT_DATABASE_FP"]
	gene_bed_fp = config["GENE_BED_FP"]
	gene_database_fp = config["GENE_DATABASE_FP"]
	eqtl_data_dir = config["EQTL_DATA_DIR"]
	expression_table_fp = config["EXPRESSION_TABLE_FP"]
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)
	
	snps = process_inputs(args.inputs,snp_database_fp,args.output_dir,args.suppress_intermediate_files)
	interactions = find_interactions(fragment_database_fp,hic_data_dir,args.include_cell_lines,args.exclude_cell_lines,args.output_dir,args.suppress_intermediate_files)
	genes = find_genes(snps,interactions,fragment_database_fp,gene_bed_fp,args.output_dir,args.suppress_intermediate_files)
	eqtls,num_sig = find_eqtls(snps,genes,eqtl_data_dir,gene_database_fp,args.fdr_threshold,args.local_databases_only,args.num_processes,args.output_dir,args.suppress_intermediate_files)
	gene_exp = get_gene_expression_info(eqtls,expression_table_fp,args.output_dir,args.suppress_intermediate_files)
	pathways = interrogate_pathways(eqtls,args.fdr_threshold,args.output_dir)
	produce_output(snps,genes,eqtls,num_sig,gene_exp,args.output_dir,args.suppress_graphs)

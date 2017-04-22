#!/usr/bin/env python

from itertools import cycle
from sets import Set
from wikipathways_api_client import WikipathwaysApiClient
import argparse,ast,bisect,configparser,csv,json,multiprocessing,os,pandas,pybedtools,re,requests,shutil,sqlite3,time

from Bio import SeqIO
from Bio import Restriction
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import Restriction_Dictionary

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.ticker import FuncFormatter

def process_inputs(inputs,snp_database_fp,fragment_database_fp,output_dir,suppress_intermediate_files=False):
	print "Processing input..."
	snp_db = sqlite3.connect(snp_database_fp)
	snp_db.text_factory = str
	snp_index = snp_db.cursor()
	fragment_index_db = sqlite3.connect(fragment_database_fp)
	fragment_index_db.text_factory = str
	fragment_index = fragment_index_db.cursor()
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
						#Query fragmentIndex to find out to which fragment the SNP belongs
						fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",[snp[1],snp[2],snp[2]])
						snp_fragment_result = fragment_index.fetchone()
						if snp_fragment_result == None:
							print "Warning: error retrieving SNP fragment for SNP " + snp[0]
						else:
							snps[snp[0]]={ "chr": snp[1], "locus": snp[2], "frag": snp_fragment_result[0] }
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
				#Query fragmentIndex to find out to which fragment the SNP belongs
				fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",[snp[1],snp[2],snp[2]])
				snp_fragment_result = fragment_index.fetchone()
				if snp_fragment_result == None:
					print "Warning: error retrieving SNP fragment for SNP " + snp[0]
				else:
					snps[snp[0]]={ "chr": snp[1], "locus": snp[2], "frag": snp_fragment_result[0] }
		
	if not suppress_intermediate_files:
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		with open(output_dir + "/snps.txt",'w') as snpfile:
			for snp,props in snps.items():
				snpfile.write("%s\t%s\t%s\t%s\n" % (snp,props["chr"],props["locus"],props["frag"]))
	return snps

def find_interactions(snps,hic_data_dir,include,exclude,output_dir,suppress_intermediate_files=False):
	print "Finding interactions..."
	#Look for all interactions involving SNP fragments in the HiC databases
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
			print "\tSearching cell line " + cell_line
			for replicate in os.listdir(hic_data_dir + '/' + cell_line):
				if replicate.endswith(".db"):
					rep_db = sqlite3.connect(hic_data_dir + '/' + cell_line + '/' + replicate)
					rep_db.text_factory = str
				        rep_ints = rep_db.cursor()
					print "\t\tSearching replicate " + replicate
					for snp in snps.keys():
						interactions[snp][cell_line] = Set([])
						print "\t\t\tFinding interactions for " + snp
						for interaction in rep_ints.execute("SELECT chr2, fragment2 FROM interactions WHERE chr1=? AND fragment1=?", [snps[snp]["chr"],snps[snp]["frag"]]):
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
				fragment_index.execute("SELECT start, end FROM fragments WHERE chr=? and fragment=?",[interaction[0],interaction[1]])
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
	num_sig = {} #Number of eQTLs deemed significant under the given threshold
	num_tests = 0 #Total number of tests done
        num_tests += query_local_databases(eqtl_data_dir,genes,eqtls,p_values)
	if not local_databases_only:
		num_tests += query_GTEx_service(snps,genes,eqtls,p_values,num_processes,output_dir)
	process_eqtls(snps,genes,eqtls,gene_database_fp)
	eqtlfile = None
	if not suppress_intermediate_files:
		eqtlfile = open(output_dir + '/eqtls.txt','w')
	for snp in eqtls.keys():
		num_sig[snp] = 0
		for gene in eqtls[snp].keys():
			if not gene == "snp_info":
				for tissue in eqtls[snp][gene]["tissues"].keys():
					eqtls[snp][gene]["tissues"][tissue]["qvalue"] = compute_fdr(eqtls[snp][gene]["tissues"][tissue]["pvalue"],p_values,num_tests)
					if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
						num_sig[snp] += 1
					if not suppress_intermediate_files:
						eqtlfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snp,eqtls[snp]["snp_info"]["chr"],eqtls[snp]["snp_info"]["locus"],gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],eqtls[snp][gene]["cell_lines"],eqtls[snp][gene]["cis?"],eqtls[snp][gene]["p_thresh"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"],eqtls[snp][gene]["tissues"][tissue]["qvalue"],eqtls[snp][gene]["tissues"][tissue]["effect_size"]))
	return (eqtls,num_sig)

def query_local_databases(eqtl_data_dir,genes,eqtls,p_values):
	print "\tQuerying local databases."
	num_tests = 0
	for db in os.listdir(eqtl_data_dir): #Iterate through databases of eQTLs by tissue type
		tissue = db[:db.rfind('.')]
		print "\t\tQuerying " + tissue
		eqtl_index_db = sqlite3.connect(eqtl_data_dir + '/' + db)
		eqtl_index_db.text_factory = str
		eqtl_index = eqtl_index_db.cursor()
		for snp in genes.keys():
			for gene in genes[snp].keys():
				num_tests += 1
				for eqtl in eqtl_index.execute("SELECT pvalue, effect_size FROM eqtls WHERE rsID=? AND gene_name=?",(snp,gene)): #Pull down all eQTLs related to a given SNP to test for relevance:
					if not eqtls.has_key(snp):
						eqtls[snp] = {}
					if not eqtls[snp].has_key(gene):
						eqtls[snp][gene] = {"tissues": {}}
					eqtls[snp][gene]["tissues"][tissue] = {"pvalue": eqtl[0], "effect_size": eqtl[1] }
					bisect.insort(p_values,eqtl[0])
	return num_tests

def query_GTEx_service(snps,genes,eqtls,p_values,num_processes,output_dir):
	tissues = Set(["Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Transformed_fibroblasts","Colon_Sigmoid","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood"])
	if num_processes > 10:
		num_processes = 10
	manager = multiprocessing.Manager()
	num_tests = 0
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
	print "\tQuerying GTEx online service..."
	gtexResponses = manager.list()
	procPool = multiprocessing.Pool(processes=num_processes)
	for i,reqList in enumerate(reqLists,start=1):
		procPool.apply_async(send_GTEx_query, (i,len(reqLists),reqList,gtexResponses))
		time.sleep(10)
	procPool.close()
	procPool.join()
	print "\t\tNumber of GTEx responses received: " + str(len(gtexResponses))
	results = []
	failed_requests = []
	for response in gtexResponses:
		try:
			results += response[1].json()["result"]
		except Exception as e:
			print "\t\tWARNING: bad response (%s)" % response[1]
			failed_requests.append(response[0])
	if failed_requests:
		with open(output_dir + "/failed_GTEx_requests.txt",'w') as failed_requests_file:
			failed_requests_file.write(str(failed_requests) + '\n')
	for result in results:
		if (str(result["geneSymbol"]) == "gene not found" or not snps.has_key(result["snpId"])) or result["pvalue"] == "NA":
			continue
		num_tests += 1
		snp = result["snpId"]
		gene = result["geneSymbol"]
		if not eqtls.has_key(snp):
			eqtls[snp] = {}
		if not eqtls[snp].has_key(gene):
			eqtls[snp][gene] = {"tissues": {}}
		p = float(result["pvalue"])
		effect_size = float(result["beta"])
		eqtls[snp][gene]["tissues"][result["tissueId"]] = { "pvalue": p, "effect_size": effect_size }
		bisect.insort(p_values,p)
	return num_tests

def process_eqtls(snps,genes,eqtls,gene_database_fp):
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	col_names = [res[1] for res in gene_index.execute("PRAGMA table_info(genes)")]
	include_p_thresh = "p_thresh" in col_names
	if include_p_thresh:
		query_str = "SELECT chr, start, end, p_thresh FROM genes WHERE symbol=?"
	else:
		query_str = "SELECT chr, start, end FROM genes WHERE symbol=?"
	for snp in eqtls.keys():
		eqtls[snp]["snp_info"] = { "chr": snps[snp]["chr"], "locus": snps[snp]["locus"] }
		for gene in eqtls[snp].keys():
				if not gene == "snp_info":
					gene_chr = "NA"
					gene_start = "NA"
					gene_end = "NA"
					max_length = 0
					for gene_stat in gene_index.execute(query_str, [gene]):
						if abs(gene_stat[2] - gene_stat[1]) > max_length: #Consider "canonical" to be the longest record where multiple records are present
							if gene_stat[0].startswith("chr"):
								gene_chr = gene_stat[0][gene_stat[0].find("chr")+3:]
							else:
								gene_chr = gene_stat[0]
							gene_start = gene_stat[1]
							gene_end = gene_stat[2]
							if include_p_thresh:
								p_thresh = gene_stat[3]
							else:
								p_thresh = "NA"
							max_length = gene_stat[2] - gene_stat[1]
					cis = gene_chr == eqtls[snp]["snp_info"]["chr"] and (eqtls[snp]["snp_info"]["locus"] > gene_start - 1000000 and eqtls[snp]["snp_info"]["locus"] < gene_end + 1000000) #eQTL is cis if the SNP is within 1Mbp of the gene
					
					eqtls[snp][gene]["gene_chr"] = gene_chr
					eqtls[snp][gene]["gene_start"] = gene_start
					eqtls[snp][gene]["gene_end"] = gene_end
					try:
						eqtls[snp][gene]["cell_lines"] = list(genes[snp][gene])
					except KeyError:
						eqtls[snp][gene]["cell_lines"] = "NA"
					eqtls[snp][gene]["p_thresh"] = p_thresh
					if cis:
						eqtls[snp][gene]["cis?"] = True
					else:
						eqtls[snp][gene]["cis?"] = False

def send_GTEx_query(num,num_reqLists,reqList,gtexResponses):
	try:
		while True:
			print "\t\tSending request %s of %s" % (num,num_reqLists)
			res = requests.post("http://gtexportal.org/api/v6p/dyneqtl?v=clversion", json=reqList)
			if res.status_code == 200:
				gtexResponses.append((reqList,res))
				time.sleep(30)
				return
			elif res.status_code == 500:
				print "\t\tThere was an error processing request %s. Writing to failed request log and continuing." % num
				gtexResponses.append((reqList,"Processing error"))
				time.sleep(30)
				return
			else:
			    print "\t\tRequest %s received response with status %s. Trying again in five minutes." % (num,res.status_code)
			    time.sleep(300)
	except requests.exceptions.ConnectionError:
		try:
			print "\t\tWarning: Request %s experienced a connection error. Retrying in five minutes." % num
			time.sleep(300)
			while True:
				print "\t\tSending request %s of %s" % (num,num_reqLists)
				res = requests.post("http://gtexportal.org/api/v6p/dyneqtl?v=clversion", json=reqList)
				print "\t\tRequest %s response: %s" % (num,res.status_code)
				if res.status_code == 200:
					gtexResponses.append((reqList,res))
					time.sleep(30)
					return
				elif res.status_code == 500:
					print "\t\tThere was an error processing request %s. Writing to failed request log and continuing." % num
					gtexResponses.append((reqList,"Processing error"))
					time.sleep(30)
					return
				else:
					print "\t\tRequest %s received response with status: %s. Trying again in five minutes." % (num,res.status_code)
					time.sleep(300)
		except requests.exceptions.ConnectionError:
			print "\t\tRetry failed. Continuing, but results will be incomplete."
			gtexResponses.append((reqList,"Connection failure"))
			time.sleep(300)
			return

def get_gene_expression_information(eqtls,expression_table_fp,output_dir):
	print "Getting gene expression information..."
	gene_df = pandas.read_table(expression_table_fp,index_col='Symbol')
	gene_exp = pandas.DataFrame(data=None,columns=gene_df.columns)
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			if not gene in gene_exp.index:
				try:
					gene_exp = gene_exp.append(gene_df.ix[gene])
				except KeyError:
					print "Warning: No gene expression information for %s" % gene
	gene_exp.to_csv(path_or_buf=output_dir+"/gene_expression_table.txt",sep='\t')

def compute_fdr(p,p_values,num_tests):
	positives = bisect.bisect_right(p_values, p)
	return float(num_tests) * p / float(positives)

def produce_summary(eqtls,expression_table_fp,output_dir):
	#TODO: Tidy up this method using pandas
	print "Producing output..."
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)
	summary = open(output_dir + "/summary.txt",'w')
	summary.write("SNP\tSNP_Chromosome\tSNP_Locus\tGene_Name\tGene_Chromosome\tGene_Start\tGene_End\tTissue\tp-value\tq-value\tEffect_Size\tCell_Lines\tGTEx_cis_p_Threshold\tcis_SNP-gene_interaction\tSNP-gene_Distance\tExpression_Level_In_eQTL_Tissue\tMax_Expressed_Tissue\tMaximum_Expression_Level\tMin_Expressed_Tissue\tMin_Expression_Level\n")
	gene_df = pandas.read_table(expression_table_fp,index_col="Description")
	gene_exp = {} #Cache of already-accessed gene expression data
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			if not gene == "snp_info":
				distance_from_snp = 0
				if(eqtls[snp][gene]["gene_chr"] == "NA"):
					distance_from_snp = "NA" #If gene location information is missing
				elif(not eqtls[snp]["snp_info"]["chr"] == eqtls[snp][gene]["gene_chr"]):
					distance_from_snp = "NA" #Not applicable to trans interactions
				elif(eqtls[snp]["snp_info"]["locus"] < eqtls[snp][gene]["gene_start"]):
					distance_from_snp = eqtls[snp][gene]["gene_start"] - eqtls[snp]["snp_info"]["locus"]
				elif(eqtls[snp]["snp_info"]["locus"] > eqtls[snp][gene]["gene_end"]):
					distance_from_snp = eqtls[snp]["snp_info"]["locus"] - eqtls[snp][gene]["gene_end"]

				for tissue in eqtls[snp][gene]["tissues"].keys():
					if not gene_exp.has_key(gene):
						gene_exp[gene] = {}
						try:
							gene_exp[gene][tissue] = gene_df.at[gene,tissue].max() #Allow for the fact that there may be more than one entry in the gene expression table for each gene
							gene_exp[gene]["max"] = gene_df.ix[gene].idxmax()
							if not isinstance(gene_exp[gene]["max"],str):
								gene_exp[gene]["max"] = gene_df.ix[gene].max().idxmax()
							gene_exp[gene][gene_exp[gene]["max"]] = gene_df.at[gene,gene_exp[gene]["max"]].max()
							gene_exp[gene]["min"] = gene_df.ix[gene].idxmin()
							if not isinstance(gene_exp[gene]["min"],str):
								gene_exp[gene]["min"] = gene_df.ix[gene].min().idxmin()
							gene_exp[gene][gene_exp[gene]["min"]] = gene_df.at[gene,gene_exp[gene]["min"]].max()
						except KeyError:
							print "\t\tWarning: No expression information for %s in %s" % (gene,tissue)
							gene_exp[gene]["max"] = "NA"
							gene_exp[gene]["min"] = "NA"
							gene_exp[gene][tissue] = "NA"
					if not gene_exp[gene].has_key(tissue):
						try:
							gene_exp[gene][tissue] = gene_df.at[gene,tissue].max()
						except KeyError:
							print "\t\tWarning: No expression information for %s in %s" % (gene,tissue)
							gene_exp[gene][tissue] = "NA"

					summary.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (snp,eqtls[snp]["snp_info"]["chr"],eqtls[snp]["snp_info"]["locus"],gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"],eqtls[snp][gene]["tissues"][tissue]["qvalue"],eqtls[snp][gene]["tissues"][tissue]["effect_size"]))
					if eqtls[snp][gene]["cell_lines"] == "NA":
						summary.write("NA")
					else:
						for i,cell_line in enumerate(eqtls[snp][gene]["cell_lines"],start=1):
							if not i == len(eqtls[snp][gene]["cell_lines"]):
								summary.write("%s," % cell_line)
							else:
								summary.write(cell_line)
					summary.write("\t%s\t%s\t%s\t%s\t%s\t" % (eqtls[snp][gene]["p_thresh"],eqtls[snp][gene]["cis?"],distance_from_snp,gene_exp[gene][tissue],gene_exp[gene]["max"]))
					if gene_exp[gene]["max"] == "NA":
						summary.write("NA\t")
					else:
						summary.write("%s\t" % gene_exp[gene][gene_exp[gene]["max"]])
					summary.write("%s\t" % gene_exp[gene]["min"])
					if gene_exp[gene]["min"] == "NA":
						summary.write("NA\n")
					else:
						summary.write("%s\n" % gene_exp[gene][gene_exp[gene]["min"]])
	summary.close()

def produce_overview(genes,eqtls,num_sig,output_dir):
	print "Producing overview..."
	stat_table = open(output_dir + "/overview.txt",'w')
	stat_table.write("SNP\tChromosome\tLocus\tTotal_SNP-gene_Pairs\tTotal_eQTLs\n")
	for snp in eqtls.keys():
		stat_table.write(snp + '\t' + eqtls[snp]["snp_info"]["chr"] + '\t' + str(eqtls[snp]["snp_info"]["locus"]) + '\t' + str(len(genes[snp])) + '\t' + str(len(eqtls[snp])-1) + '\n')
	stat_table.close()
	print "\tProducing graphs"
	style.use("ggplot")
	if not os.path.isdir(output_dir + "/plots"):
		os.mkdir(output_dir + "/plots")
	int_colours = "rgb"
	eqtl_colours = "myc"
	int_colours = cycle(int_colours)
	eqtl_colours = cycle(eqtl_colours)
	snps_by_chr = {}
	for snp in eqtls.keys():
		chrm = eqtls[snp]["snp_info"]["chr"]
		if not snps_by_chr.has_key(chrm):
			snps_by_chr[chrm] = []
		snps_by_chr[chrm].append((eqtls[snp]["snp_info"]["locus"],snp))
	
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

def retrieve_pathways(eqtls,fdr_threshold,num_processes,output_dir):
	print "Retrieving pathway information..."
	manager = multiprocessing.Manager()
	procPool = multiprocessing.Pool(processes=num_processes)
	pathways = {} #Keeps track of all pathways covered by genes with a statistically signficant eQTL relationship with a query SNP
	pwresults = manager.list()
	
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			if not gene == "snp_info":
				for tissue in eqtls[snp][gene]["tissues"].keys():
					if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
						procPool.apply_async(get_wikipathways_response, [snp,gene,tissue,pwresults])
	procPool.close()
	procPool.join()

	for response in pwresults:
		snp = response[0]
		gene = response[1]
		tissue = response[2]
		for pwres in response[3]:
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

def get_wikipathways_response(snp,gene,tissue,pwresults):
	wikipathways = WikipathwaysApiClient()
	kwargs = {
		'query': gene,
		'organism': 'http://identifiers.org/taxonomy/9606' #Homo sapiens
	}
	res = wikipathways.find_pathways_by_text(**kwargs)
	if res:
		pwresults.append([snp,gene,tissue,res])

def parse_snps_files(snps_files):
	snps = {}
	for snp_file in snps_files:
		with open(snp_file,'r') as snpfile:
			for line in snpfile:
				snp = line.strip().split('\t')
				snps[snp[0]] = { "chr": snp[1], "locus": int(snp[2]), "frag": snp[3] }
	return snps

def parse_interactions_files(interactions_files):
	interactions = {}
	for interactions_file in interactions_files:
		with open(interactions_file,'r') as intfile:
			for line in intfile:
				interaction = line.strip().split('\t')
				snp = interaction[0]
				cell_line = interaction[1]
				if not interactions.has_key(snp):
					interactions[snp] = {}
				if not interactions[snp].has_key(cell_line):
					interactions[snp][cell_line] = Set([])
				interactions[snp][cell_line].add((interaction[2],int(interaction[3])))
	return interactions

def parse_genes_files(genes_files):
	genes = {}
	for gene_file in genes_files:
		with open(gene_file,'r') as genefile:
			for line in genefile:
				gene = line.strip().split('\t')
				snp = gene[0]
				if not genes.has_key(snp):
					genes[snp] = {}
				if not genes[snp].has_key(gene[1]):
					genes[snp][gene[1]] = Set([])
				genes[snp][gene[1]].add(gene[2])
	return genes

def parse_eqtls_files(eqtls_files,fdr_threshold=None):
	print "Parsing eQTLs..."
	eqtls = {}
	p_values = [] #A sorted list of all p-values for use computing FDR
	num_tests = 0 #Total number of tests done
	num_sig = {} #Number of eQTLs deemed significant under the given threshold
	for eqtls_file in eqtls_files:
		lines = 0
		with open(eqtls_file,'r') as eqtlfile:	
			#Do line count for progress meter
			print "\tDetermining table size..."
			for i in eqtlfile:
				lines += 1
		lines = lines//100*100 #Get an approximation
		do_linecount = not lines == 0
		with open(eqtls_file,'r') as eqtlfile:
			for i,line in enumerate(eqtlfile):
				if do_linecount:
					if i % (lines/100) == 0:
						print "\t\tProcessed %d%%..." % ((float(i)/float(lines))*100)
				eqtl = line.strip().split('\t')
				snp = eqtl[0]
				gene = eqtl[3]
				gene_chr = eqtl[4]
				try:
					gene_start = int(eqtl[5])
					gene_end = int(eqtl[6])
				except ValueError:
					gene_start = eqtl[5]
					gene_end = eqtl[6]
					try:
						cell_lines = ast.literal_eval(eqtl[7])
					except ValueError:
						cell_lines = eqtl[7]
				cis = ast.literal_eval(eqtl[8])
				try:
					p_thresh = float(eqtl[9])
				except ValueError:
					p_thresh = eqtl[9]
				tissue = eqtl[10]
				p = float(eqtl[11])
				q = float(eqtl[12])
				effect_size = float(eqtl[13])
				if not eqtls.has_key(snp):
					eqtls[snp] = {}
					num_sig[snp] = 0
					eqtls[snp]["snp_info"] = { "chr": eqtl[1], "locus": int(eqtl[2]) }
				if not eqtls[snp].has_key(gene):
					eqtls[snp][gene] = {}
					eqtls[snp][gene]["gene_chr"] = gene_chr
					eqtls[snp][gene]["gene_start"] = gene_start
					eqtls[snp][gene]["gene_end"] = gene_end
					eqtls[snp][gene]["cell_lines"] = cell_lines
					eqtls[snp][gene]["p_thresh"] = p_thresh
					eqtls[snp][gene]["tissues"] = {}
					eqtls[snp][gene]["cis?"] = cis
					eqtls[snp][gene]["tissues"][tissue] = {"pvalue": p, "effect_size": effect_size}
				if fdr_threshold:
					bisect.insort(p_values,p)
					num_tests += 1
				else:
					eqtls[snp][gene]["tissues"][tissue]["qvalue"] = q
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			if not gene == "snp_info":
				for tissue in eqtls[snp][gene]["tissues"].keys():
					if fdr_threshold:
						eqtls[snp][gene]["tissues"][tissue]["qvalue"] = compute_fdr(eqtls[snp][gene]["tissues"][tissue]["pvalue"],p_values,num_tests)
						if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
							num_sig[snp] += 1
	return (eqtls,num_sig)

def build_snp_index(snp_dir,output_fp,config,id_col=4,chr_col=1,locus_col=2,do_not_tidy_up=False):
	if not output_fp:
		output_fp = config["SNP_DATABASE_FP"]

	print "Building SNP index..."
	if not os.path.isdir(snp_dir):
		print "Error: argument to build SNP index must be a directory."
		return

	if os.path.isfile(output_fp):
		upsert = raw_input("WARNING: Upserting input to existing SNP database (%s). Continue? [y/N] " % output_fp)
		if not upsert.lower() == 'y':
			print "Did not write to existing SNP database."
			return

	if not os.path.isdir(os.path.dirname(output_fp)):
		os.makedirs(os.path.dirname(output_fp))

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
	if not do_not_tidy_up:
		print "Tidying up..."
		shutil.rmtree(snp_dir)

def build_hic_index(input_hic_fp,output_fp=None,chr1_col=3,chr2_col=7,frag1_col=5,frag2_col=9,mapq1_col=10,mapq2_col=11,mapq_cutoff=150,do_not_tidy_up=False):
	if not output_fp:
		if not input_hic_fp.rfind('.') == -1:
			output_fp = input_hic_fp[:input_hic_fp.rfind('.')] + ".db"
		else:
			output_fp = input_hic_fp + ".db"

	if os.path.isfile(output_fp):
		overwrite = raw_input("WARNING: Overwriting existing HiC database (%s). Continue? [y/N] " % output_fp)
		if not upsert.lower() == 'y':
			print "Did not overwrite existing HiC database."
			return
		os.remove(output_fp)

	#Do line count for progress meter
	print "Determining table size..."
	hic_table = open(input_hic_fp,'r')
	lines = 0
	for i in hic_table:
	    lines += 1
	hic_table.close()
	lines = lines//100*100 #Get an approximation
	do_linecount = not lines == 0
	
	int_db = sqlite3.connect(output_fp)
	interactions = int_db.cursor()
	interactions.execute("CREATE TABLE IF NOT EXISTS interactions (chr1 text, fragment1 text, chr2 text, fragment2 text)")
	interactions.execute("CREATE INDEX IF NOT EXISTS i_1 ON interactions (chr1, fragment1)")
	chr1_col = chr1_col - 1
	chr2_col = chr2_col - 1
	frag1_col = frag1_col - 1
	frag2_col = frag2_col - 1
	mapq1_col = mapq1_col - 1
	mapq2_col = mapq2_col - 1
	with open(input_hic_fp,'r') as rao_table:
		print "Indexing HiC interaction table..."
		for i,line in enumerate(rao_table):
			if do_linecount:
				if i % (lines/100) == 0:
					print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
			interaction = line.strip().split(' ')
			if int(interaction[mapq1_col]) >= mapq_cutoff and int(interaction[mapq2_col]) >= mapq_cutoff:
				chr1 = interaction[chr1_col]
				frag1 = interaction[frag1_col]
				chr2 = interaction[chr2_col]
				frag2 = interaction[frag2_col]
				interactions.execute("INSERT INTO interactions VALUES(?,?,?,?)", [chr1,frag1,chr2,frag2])
				interactions.execute("INSERT INTO interactions VALUES(?,?,?,?)", [chr2,frag2,chr1,frag1])
        int_db.commit()
        interactions.close()
	print "Done indexing HiC interaction table."
	if not do_not_tidy_up:
		print "Tidying up..."
		os.remove(input_hic_fp)

def digest_genome(genome,restriction_enzyme,output_fp,output_db,do_not_index=False,linear=False):
	if not output_fp:
		if not genome.rfind('.') == -1:
			output_fp = "%s.%s.fragments.bed" % (genome[:genome.rfind('.')],restriction_enzyme)
		else:
			output_fp = "%s.%s.fragments.bed" % (genome,restriction_enzyme)
	if not os.path.isdir(os.path.dirname(output_fp)):
		os.makedirs(os.path.dirname(output_fp))
	if os.path.isfile(output_fp):
		overwrite = raw_input("WARNING: Overwriting existing fragment BED %s. Continue? [y/N] " % output_fp)
		if not overwrite.lower() == 'y':
			print "Did not overwrite existing fragment BED."
			return
		os.remove(output_fp)

	print "Digesting"
	if "fasta" in genome or "fa" in genome:
		genome = SeqIO.parse(open(genome,"rU"), format='fasta')
	else:
		genome = SeqIO.parse(open(genome,"rU"), format='genbank')

	with open(output_fp, 'w') as bedfile:
		output = csv.writer(bedfile,delimiter="\t")
		for chromosome in genome:
			print chromosome.id, chromosome.name
			#Digest the sequence data and return the cut points
			fragment_num = 0
			enzyme = RestrictionBatch([restriction_enzyme])
			for enzyme, cutpoints in enzyme.search(chromosome.seq, linear=linear).items():
				cutpoints.insert(0,0)
				#Covert cut points to fragments
				for index, point in enumerate(cutpoints):
					#Adjust for start and end of sequence and the offset for the cutpoint
					#ATTENTION: This will only work for MspI I will have to spend some time to make it compatible with 
					#		   any restriction enzyme such as ones with blunt ends
					#startpoint = 0 if cutpoints[index] - 1 < 0 else cutpoints[index] - 1
					startpoint = 0 if cutpoints[index] - 1 < 0 else cutpoints[index]
					endpoint = len(chromosome.seq) if index+1 >= len(cutpoints) else cutpoints[index+1] - 1

					accession = chromosome.id
					version = ''
					if "." in chromosome.id:
						accession, version = chromosome.id.split(".")
					if not accession.startswith("chr"):
						accession = "chr" + accession
					output.writerow([accession, startpoint, endpoint, fragment_num])
					fragment_num += 1
	if not do_not_index:
		build_fragment_index(output_fp,output_db)

def build_fragment_index(fragment_fp,output_db):
	if not output_db:
		if not fragment_fp.rfind('.') == -1:
			output_db = fragment_fp[:fragment_fp.rfind('.')] + ".db"
		else:
			output_db = fragment_fp + ".db"
	if os.path.isfile(output_db):
		overwrite = raw_input("WARNING: Overwriting existing fragment database %s. Continue? [y/N] " % output_db)
		if not overwrite.lower() == 'y':
			print "Did not overwrite existing fragment database."
			return
		os.remove(output_db)

	if not os.path.isdir(os.path.dirname(output_db)):
		os.path.makedirs(os.path.dirname(output_db))

	fragment_index_db = sqlite3.connect(output_db)
	fragment_index = fragment_index_db.cursor()
	fragment_index.execute("CREATE TABLE fragments (chr text, start integer, end integer, fragment integer)")
	fragment_index.execute("CREATE INDEX f_index ON fragments (chr,fragment)")
	
	with open(fragment_fp,'r') as fragments_bed:
		for line in fragments_bed:
			fragment = line.strip().split('\t')
			fragment_index.execute("INSERT INTO fragments VALUES (?,?,?,?)", [fragment[0][fragment[0].find("chr")+3:],int(fragment[1]),int(fragment[2]),fragment[3]])
	fragment_index_db.commit()

def build_gene_index(gene_files,output_bed,output_db,config,symbol_col=27,chr_col=29,start_col=30,end_col=31,p_thresh_col=None,no_header=False,do_not_tidy_up=False):
	genes = {}

	symbol_col -= 1
	chr_col -= 1
	start_col -= 1
	end_col -= 1
	if p_thresh_col:
		p_thresh_col -= 1

	if not output_bed:
		output_bed = os.path.join(config.get("Defaults","LIB_DIR"),"gene_reference.bed")

	if not output_db:
		output_bed = os.path.join(config.get("Defaults","LIB_DIR"),"gene_reference.db")

	append_bed = False
	overwrite_bed = True
	if os.path.isfile(output_bed):
		upsert = raw_input("WARNING: Appending input to existing BED file (%s). Continue? [y/N] " % output_bed)
		if not upsert.lower() == 'y':
			print "Did not append to existing gene database."
			overwrite_bed = False
		else:
			append_bed = True

	if (overwrite_bed or append_bed) and not os.path.isdir(os.path.dirname(output_bed)):
		os.makedirs(os.path.dirname(output_bed))

	upsert_db = True
	if os.path.isfile(output_db):
		upsert = raw_input("WARNING: Upserting input to existing gene database (%s). Continue? [y/N] " % output_db)
		if not upsert.lower() == 'y':
			print "Did not write to existing gene database."
			upsert_db = False

	if upsert_db and not os.path.isdir(os.path.dirname(output_db)):
		os.makedirs(os.path.dirname(output_db))

	if (not (overwrite_bed or append_bed)) and not upsert_db:
		print "No action performed; exiting."
		return

	if upsert_db:
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
	if overwrite_bed:
		bed_out = open(output_bed,'w')
	elif append_bed:
		bed_out = open(output_bed,'a')
	for gene in genes.keys():
		if bed_out:
			bed_out.write('chr%s\t%s\t%s\t%s\n' % (genes[gene]["chr"],genes[gene]["start"],genes[gene]["end"],gene))
		if upsert_db:
			if p_thresh_col:
				gene_index.execute("INSERT INTO genes VALUES (?,?,?,?,?)", [gene,genes[gene]["chr"],genes[gene]["start"],genes[gene]["end"],genes[gene]["p_thresh"]])
			else:
				gene_index.execute("INSERT INTO genes VALUES (?,?,?,?)", [gene,genes[gene]["chr"],genes[gene]["start"],genes[gene]["end"]])
	if bed_out:
		bed_out.close()
	if upsert_db:
		gene_index_db.commit()
		gene_index.close()

	if not do_not_tidy_up:
		print "Tidying up..."
		for gene_file in gene_files:
			os.remove(gene_file)

def build_eqtl_index(table_fp,output_fp=None,snp_col=23,gene_symbol_col=27,gene_chr_col=29,gene_start_col=30,gene_stop_col=31,p_val_col=6,effect_size_col=3,do_not_tidy_up=False):
	if not output_fp:
		if not table_fp.rfind('.') == -1:
			output_fp = table_fp[:table_fp.rfind('.')] + ".db"
		else:
			output_fp = table_fp + ".db"

	if not os.path.isdir(os.path.dirname(output_fp)):
		os.makedirs(os.path.dirname(output_fp))
	
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
        p_val_col -= 1
	effect_size_col -= 1
	table_index_db = sqlite3.connect(output_fp)
	table_index = table_index_db.cursor()
	table_index.execute("CREATE TABLE IF NOT EXISTS eqtls (rsID text, gene_name text, gene_chr text, gene_start integer, gene_end integer, pvalue real, effect_size real)")
	table_index.execute("CREATE INDEX IF NOT EXISTS id ON eqtls (rsID)")
	##Do line count for progress meter
	do_linecount = True
	print "Determining table size..."
        with open(table_fp,'r') as eqtl_table:
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
	if not do_not_tidy_up:
		print "Tidying up..."
		os.remove(table_fp)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',required=True,help="The the dbSNP IDs or loci of SNPs of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-c","--config",default=os.path.join(os.path.dirname(__file__),"../docs/codes3d.conf"),help="The configuration file to be used for this hiCquery run (default: conf.py)")
	parser.add_argument("-n","--include_cell_lines",nargs='+',help="Space-separated list of cell lines to include (others will be ignored). NOTE: Mutually exclusive with EXCLUDE_CELL_LINES.")
	parser.add_argument("-x","--exclude_cell_lines",nargs='+',help="Space-separated list of cell lines to exclude (others will be included). NOTE: Mutually exclusive with INCLUDE_CELL_LINES.")
	parser.add_argument("-o","--output_dir",default="codes3d_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	parser.add_argument("-l","--local_databases_only",action="store_true",default=False,help="Consider only local databases. Will only include cis-eQTLs if using downloadable GTEx dataset.")
	parser.add_argument("-s","--suppress_intermediate_files",action="store_true",default=False,help="Do not produce intermediate files. These can be used to run the pipeline from an intermediate stage in the event of interruption.")
	parser.add_argument("-p","--num_processes",type=int,default=1,help="Desired number of processes for multithreading (default: 1).")
	parser.add_argument("-f","--fdr_threshold",type=float,default=0.05,help="The FDR threshold to consider an eQTL statistically significant (default: 0.05).")
	args = parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config)
	snp_database_fp = config.get("Defaults","SNP_DATABASE_FP")
	hic_data_dir = config.get("Defaults","HIC_DATA_DIR")
	fragment_bed_fp = config.get("Defaults","FRAGMENT_BED_FP")
	fragment_database_fp = config.get("Defaults","FRAGMENT_DATABASE_FP")
	gene_bed_fp = config.get("Defaults","GENE_BED_FP")
	gene_database_fp = config.get("Defaults","GENE_DATABASE_FP")
	eqtl_data_dir = config.get("Defaults","EQTL_DATA_DIR")
	expression_table_fp = config.get("Defaults","EXPRESSION_TABLE_FP")
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)
	
	snps = process_inputs(args.inputs,snp_database_fp,fragment_database_fp,args.output_dir,args.suppress_intermediate_files)
	interactions = find_interactions(snps,hic_data_dir,args.include_cell_lines,args.exclude_cell_lines,args.output_dir,args.suppress_intermediate_files)
	genes = find_genes(interactions,fragment_database_fp,gene_bed_fp,args.output_dir,args.suppress_intermediate_files)
	eqtls,num_sig = find_eqtls(snps,genes,eqtl_data_dir,gene_database_fp,args.fdr_threshold,args.local_databases_only,args.num_processes,args.output_dir,suppress_intermediate_files=args.suppress_intermediate_files)
	produce_summary(eqtls,expression_table_fp,args.output_dir)
	produce_overview(genes,eqtls,num_sig,args.output_dir)
	pathways = retrieve_pathways(eqtls,args.fdr_threshold,args.num_processes,args.output_dir)

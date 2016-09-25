#!/usr/bin/python
from configobj import ConfigObj
from itertools import cycle
from sets import Set
from wikipathways_api_client import WikipathwaysApiClient
import argparse,ast,bisect,json,multiprocessing,os,pandas,pybedtools,re,requests,sqlite3,time

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
						fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",["chr" + snp[1],snp[2],snp[2]])
						snp_fragment_result = fragment_index.fetchone()
						if snp_fragment_result == None:
							print "Warning: error retrieving SNP fragment for SNP " + snp
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
				fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",["chr" + snp[1],snp[2],snp[2]])
				snp_fragment_result = fragment_index.fetchone()
				if snp_fragment_result == None:
					print "Warning: error retrieving SNP fragment for SNP " + snp
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
						eqtlfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snp,eqtls[snp]["snp_info"]["chr"],eqtls[snp]["snp_info"]["locus"],gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],eqtls[snp][gene]["cell_lines"],eqtls[snp][gene]["cis?"],eqtls[snp][gene]["p_thresh"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"]))
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
				for eqtl in eqtl_index.execute("SELECT pvalue FROM eqtls WHERE rsID=? AND gene_symbol=?",(snp,gene)): #Pull down all eQTLs related to a given SNP to test for relevance:
					if not eqtls.has_key(snp):
						eqtls[snp] = {}
					if not eqtls[snp].has_key(gene):
						eqtls[snp][gene] = {"tissues": {}}
					eqtls[snp][gene]["tissues"][tissue] = {"pvalue": eqtl[0] }
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
		eqtls[snp][gene]["tissues"][result["tissueId"]] = { "pvalue": p }
		bisect.insort(p_values,p)
	return num_tests

def process_eqtls(snps,genes,eqtls,gene_database_fp):
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	for snp in eqtls.keys():
		eqtls[snp]["snp_info"] = { "chr": snps[snp]["chr"], "locus": snps[snp]["locus"] }
		for gene in eqtls[snp].keys():
				if not gene == "snp_info":
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
	summary.write("SNP\tSNP_Chromosome\tSNP_Locus\tGene_Name\tGene_Chromosome\tGene_Start\tGene_End\tTissue\tp-value\tq-value\tCell_Lines\tGTEx_cis_p_Threshold\tcis_SNP-gene_interaction\tSNP-gene_Distance\tExpression_Level_In_eQTL_Tissue\tMax_Expressed_Tissue\tMaximum_Expression_Level\tMin_Expressed_Tissue\tMin_Expression_Level\n")
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

					summary.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (snp,eqtls[snp]["snp_info"]["chr"],eqtls[snp]["snp_info"]["locus"],gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"],eqtls[snp][gene]["tissues"][tissue]["qvalue"]))
					if eqtls[snp][gene]["cell_lines"] == "NA":
						summary.write("NA")
					else:
						for i,cell_line in enumerate(eqtls[snp][gene]["cell_lines"],start=1):
							if not i == len(eqtls[snp][gene]["cell_lines"]):
								summary.write(cell_line + ',')
							else:
								summary.write(cell_line)
					summary.write("\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (eqtls[snp][gene]["p_thresh"],eqtls[snp][gene]["cis?"],distance_from_snp,gene_exp[gene][tissue],gene_exp[gene]["max"],gene_exp[gene][gene_exp[gene]["max"]],gene_exp[gene]["min"],gene_exp[gene][gene_exp[gene]["min"]]))
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
	for interaction_file in interactions_files:
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

def parse_eqtls_files(eqtls_files,fdr_threshold):
	eqtls = {}
	p_values = [] #A sorted list of all p-values for use computing FDR
	num_tests = 0 #Total number of tests done
	num_sig = {} #Number of eQTLs deemed significant under the given threshold
	for eqtls_file in eqtls_files:
		with open(eqtls_file,'r') as eqtlfile:
			for line in eqtlfile:
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
				eqtls[snp][gene]["tissues"][tissue] = {"pvalue": p}
				bisect.insort(p_values,p)
				num_tests += 1
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			if not gene == "snp_info":
				for tissue in eqtls[snp][gene]["tissues"].keys():
					eqtls[snp][gene]["tissues"][tissue]["qvalue"] = compute_fdr(eqtls[snp][gene]["tissues"][tissue]["pvalue"],p_values,num_tests)
					if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
						num_sig[snp] += 1
	return (eqtls,num_sig)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',required=True,help="The the dbSNP IDs or loci of SNPs of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-c","--config",default="docs/conf.py",help="The configuration file to be used for this hiCquery run (default: conf.py)")
	parser.add_argument("-n","--include_cell_lines",nargs='+',help="Space-separated list of cell lines to include (others will be ignored). NOTE: Mutually exclusive with EXCLUDE_CELL_LINES.")
	parser.add_argument("-x","--exclude_cell_lines",nargs='+',help="Space-separated list of cell lines to exclude (others will be included). NOTE: Mutually exclusive with INCLUDE_CELL_LINES.")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	parser.add_argument("-l","--local_databases_only",action="store_true",default=False,help="Consider only local databases. Will only include cis-eQTLs if using downloadable GTEx dataset.")
	parser.add_argument("-s","--suppress_intermediate_files",action="store_true",default=False,help="Do not produce intermediate files. These can be used to run the pipeline from an intermediate stage in the event of interruption.")
	parser.add_argument("-p","--num_processes",type=int,default=1,help="Desired number of processes for multithreading (default: 1).")
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
	
	snps = process_inputs(args.inputs,snp_database_fp,fragment_database_fp,args.output_dir,args.suppress_intermediate_files)
	interactions = find_interactions(snps,hic_data_dir,args.include_cell_lines,args.exclude_cell_lines,args.output_dir,args.suppress_intermediate_files)
	genes = find_genes(interactions,fragment_database_fp,gene_bed_fp,args.output_dir,args.suppress_intermediate_files)
	eqtls,num_sig = find_eqtls(snps,genes,eqtl_data_dir,gene_database_fp,args.fdr_threshold,args.local_databases_only,args.num_processes,args.output_dir,suppress_intermediate_files=args.suppress_intermediate_files)
	produce_summary(eqtls,expression_table_fp,args.output_dir)
	produce_overview(genes,eqtls,num_sig,args.output_dir)
	pathways = retrieve_pathways(eqtls,args.fdr_threshold,args.num_processes,args.output_dir)

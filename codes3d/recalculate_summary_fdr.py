#!/usr/bin/python

import os
import csv
import argparse
import ast
import codes3d
import configparser

def get_eqtls(summary_file, eqtl_file):
    fdr_threshold = None
    interactions = []
    summary = []
    eqtl_list = []
    eqtls = {}
    p_values = [] #A sorted list of all p-values for use computing FDR
    num_tests = 0 #Total number of tests done
    num_sig = {} #Number of eQTLs deemed significant under the given threshold
    summ_file = open(summary_file, 'r')
    summ_reader = csv.reader(summ_file, delimiter = '\t')
    for line in summ_reader:
        if not line[0].startswith('SNP'):
            pair = (line[0], line[3])
            if pair not in interactions:
                interactions.append(pair)
            summary.append(line)
    eqtl_file = open(eqtl_file, 'r')
    eqtl_reader = csv.reader(eqtl_file, delimiter = '\t')
    for line in eqtl_reader:
        pair = (line[0], line[3])
        if pair in interactions:
            eqtl_list.append(line)
    lines = len(eqtl_list)//100*100
    do_linecount = not lines == 0
    for i, line in enumerate(eqtl_list):
        if do_linecount:
            if i % (lines/100) == 0:
                pass
                #print('\t\tProcessed %d%%...' %((float(i)/lines)*100))
        eqtl = line
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
        try:
            q = float(eqtl[12])
            effect_size = float(eqtl[13])
        except IndexError:
            q = 0.0 #@Tayaza: Handle previous versions without qvalue 
            effect_size = 0.0 #@Tayaza: Handle previous versions 
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
        eqtls[snp][gene]["tissues"][tissue] = {"pvalue": p, "effect_size": effect_size} \
            #@Tayaza: To capture all eQTL tissues and not just the first one
        if fdr_threshold:
            bisect.insort(p_values,p)
            num_tests += 1
        else:
            eqtls[snp][gene]["tissues"][tissue]["qvalue"] = q
        
    i = 0
    for snp in eqtls.keys():
        for gene in eqtls[snp].keys():
            if not gene == "snp_info":
                for tissue in eqtls[snp][gene]["tissues"].keys():
                    if fdr_threshold:
                        eqtls[snp][gene]["tissues"][tissue]["qvalue"] \
                            = compute_fdr(eqtls[snp][gene]["tissues"][tissue]["pvalue"],p_values,num_tests)
                        if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
                            num_sig[snp] += 1
                    qvalue = eqtls[snp][gene]["tissues"][tissue]["qvalue"]
                    pvalue = eqtls[snp][gene]["tissues"][tissue]["pvalue"]
                    if qvalue <= 0.05:
                        i += 1
    return(eqtls)

def compute_fdr(p,p_values,num_tests):
    positives = bisect.bisect_right(p_values, p)
    return float(num_tests) * p / float(positives)
        
        

if __name__=='__main__':
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument('-s', '--summary_file', required=True, help='CoDeS3D output summary.txt')
    parser.add_argument('-e', '--eqtl_file', required=True, help='CoDeS3D eqtls.txt')
    parser.add_argument("-o","--output_dir",default="codes3d_summary",\
                            help="The directory in which to output results (\"codes3d_output\" by default).")
    parser.add_argument("-c","--config",default="docs/codes3d.conf",\
                            help="The configuration file to be used in this instance (default: conf.py)")
    parser.add_argument("-f","--fdr_threshold",type=float,default=0.05,\
                            help="The FDR threshold to consider an eQTL statistically significant (default: 0.05).")
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    expression_table_fp = config.get("Defaults","EXPRESSION_TABLE_FP")
    eqtls = get_eqtls(args.summary_file, args.eqtl_file)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    print(expression_table_fp)
    codes3d.produce_summary(eqtls,expression_table_fp,args.output_dir)

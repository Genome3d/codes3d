#usr/bin/env python

from itertools import cycle
import wikipathways_api_client
import argparse
import ast
import bisect
import configparser
import csv
import json
import multiprocessing
import operator
import os
import sys
import pandas
import psutil
import pybedtools
import re
import requests
import shutil
import sqlite3
import time
from operator import itemgetter
import Bio
from Bio import SeqIO
from Bio import Restriction
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import Restriction_Dictionary
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.ticker import FuncFormatter
import rpy2.robjects as R

def parse_parameters(restriction_enzymes, include_cell_line, exclude_cell_line):
    """Validate user parameters -r, -n and -x.

    Args:
        restriction_enzymes: space-delimited list of restriction enzymes from
            user. Limits program to Hic libraries restricted by specified enzyme.
        include_cell_line: space-delimited list of cell_lines from -n.
        exclude_cell_line: space-delimited list of cell_lines from -x
    Returns:
        res_enzymes: A list of restriction enzymes of which HiC libraries are
            included.
        include_cells: A list of validated HiC cell lines to be querried.
        exclude_cells: A list of validated HiC cell lines to be excluded.
    """
    print('Parsing HiC library restriction enzymes...')
    res_enzymes = []
    include_cells = []
    exclude_cells = []
    if restriction_enzymes == None:
        res_enzymes = HIC_RESTRICTION_ENZYMES
    else:
        for enzyme in restriction_enzymes:
            # Handle user input case issues
            enzymes_upper = [e.upper() for e in HIC_RESTRICTION_ENZYMES]
            if not enzyme.upper() in enzymes_upper:
                print('Warning: NO HiC libraries restricted with %s.' % enzyme)
            else:
                res_enzymes.append(HIC_RESTRICTION_ENZYMES[
                    enzymes_upper.index(enzyme.upper())])
    if not res_enzymes:
        sys.exit("Program terminating: No HiC libraries are prepared with " +
              "your restriction enzyme(s).")
    cell_lines = []
    for enzyme in res_enzymes:
        cell_lines += os.listdir(os.path.join(lib_dir, enzyme, 'hic_data'))
    cell_lines_upper = [c.upper() for c in cell_lines]
    to_help = '1.\t' + cell_lines[0]
    if len(cell_lines) > 1:
        for i, cell in enumerate(cell_lines, start=1):
            to_help += '\n' + str(i+1) + '.\t' + cell
    if include_cell_line == None:
        include_cells = None
    else:
        for cell in include_cell_line:
                if cell.upper() in cell_lines_upper:
                    include_cells.append(cell_lines[cell_lines_upper.index(cell.upper())])
        if include_cells != None and len(include_cells) == 0:
            sys.exit("We cannot find one or more of the cell lines you " +\
                     "passed to the -n option. Following is the list of " +\
                     "cell lines we have: \n" + to_help)
    if exclude_cell_line == None:
        exclude_cells = None
    else:
        for cell in exclude_cell_line:
            if cell.upper() in cell_lines_upper:
                exclude_cells.append(cell_lines[cell_lines_upper.index(cell.upper())])
        if exclude_cells != None and len(exclude_cells) == 0:
            sys.exit("We cannot find one or more of the cell lines you " +\
                     "passed to the -x option. Following is the list of " +\
                     "cell lines we have: \n" + to_help)

    return res_enzymes, include_cells, exclude_cells

def process_inputs(inputs, snp_database_fp, lib_dir,
                   restriction_enzymes, output_dir,
                   suppress_intermediate_files=False):
    """Finds restriction fragments in which SNPs lie.

    Args:
        inputs: File(s) (or stdin) containing SNP rsIDs or genomic positions
          in the format chr<x>:<locus>
        snp_database_fp: ../../lib/snp_index_dbSNP_b147.db
        lib_dir: ../lib #To point to fragment databases of different hic restrictions.
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
        A dict named snps with the ff structure:
            {'rs9462794':{
                          'chr': '6',
                          'locus': 12445846,
                          'fragments':[
                                       {'frag': 30968, 'enzyme': 'MboI'},
                                       {'frag': 3664, 'enzyme': 'HindIII'}
                                    ]
                        },
             'rs12198798':{...},
             'rs6909834':{...}
             }

        If suppress_intermediate_files=False, a snps.txt file is written
          with the ff columns:
            1. SNP rsID
            2. SNP chromosome
            3. SNP position
            4. Fragment ID
            5. Fragment restriction enzyme
    """
    print("Processing SNP input...")
    snp_db = sqlite3.connect(snp_database_fp)
    snp_db.text_factory = str
    snp_index = snp_db.cursor()
    #fragment_index_db = sqlite3.connect(fragment_database_fp)
    #fragment_index_db.text_factory = str
    #fragment_index = fragment_index_db.cursor()
    snps = {}
    for input in inputs:
        if os.path.isfile(input):
            infile = open(input, 'r')
            for line in infile:
                id = line.strip().split(' ')[0]
                snp = None
                if id.startswith('rs'):
                    snp_index.execute("SELECT * FROM snps WHERE rsID=?", (id,))
                else:
                    if not id.strip():
                        continue
                    chr = id[id.find("chr") + 3:id.find(':')]
                    locus = int(id[id.find(':') + 1:])
                    snp_index.execute(
                        "SELECT * FROM snps WHERE chr=? and locus=?", [chr, locus])
                snp = snp_index.fetchone()
                if snp is None:
                    print("Warning: %s does not exist in SNP database." % id)
                else:
                    # Query each fragment_bed_fp for SNP fragment.
                    for enzyme in restriction_enzymes:
                        fragment_database_fp = os.path.join(lib_dir,enzyme,'dna.fragments.db')
                        fragment_index_db = sqlite3.connect(fragment_database_fp)
                        fragment_index_db.text_factory = str
                        fragment_index = fragment_index_db.cursor()
                        fragment_index.execute("SELECT fragment FROM fragments " +
                                           "WHERE chr=? AND start<=? AND end>=?",
                                           [snp[1], snp[2], snp[2]])
                        snp_fragment_result = fragment_index.fetchone()
                        if snp_fragment_result is None:
                            print("Warning: error retrieving SNP fragment for SNP " +
                                  snp[0])
                        else:
                            if not snp[0] in snps:
                                snps[snp[0]] = {"chr": '', "locus": "", "fragments": []}
                                snps[snp[0]]["chr"] = snp[1]
                                snps[snp[0]]["locus"] = snp[2]
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                               "enzyme": enzyme})
                            else:
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                               "enzyme": enzyme})
            infile.close()
        else:
            snp = None
            if input.startswith("rs"):
                snp_index.execute("SELECT * FROM snps WHERE rsID=?", (input,))
            else:
                chr = input[input.find("chr") + 3:input.find(':')]
                locus = int(input[input.find(':') + 1:])
                snp_index.execute("SELECT * FROM snps WHERE chr=? and locus=?",
                                  [chr, locus])
                snp = snp_index.fetchone()
                if snp is None:
                    print(
                        "Warning: %s does not exist in SNP database." %
                        input)
                else:
                    # Query each fragment_bed_fp for SNP fragment.
                    for enzyme in restriction_enzymes:
                        fragment_database_fp = os.path.join(lib_dir,enzyme,'dna.fragments.db')
                        fragment_index_db = sqlite3.connect(fragment_database_fp)
                        fragment_index_db.text_factory = str
                        fragment_index = fragment_index_db.cursor()
                        fragment_index.execute("SELECT fragment FROM fragments " +
                                           "WHERE chr=? AND start<=? AND end>=?",
                                           [snp[1], snp[2], snp[2]])
                        snp_fragment_result = fragment_index.fetchone()
                        if snp_fragment_result is None:
                            print("Warning: error retrieving SNP fragment for SNP " +
                                  snp[0])
                        else:
                            if not snp[0] in snps:
                                snps[snp[0]] = {"chr": '', "locus": "", "fragments": []}
                                snps[snp[0]]["chr"] = snp[1]
                                snps[snp[0]]["locus"] = snp[2]
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                               "enzyme": enzyme})
                            else:
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                               "enzyme": enzyme})
    if not suppress_intermediate_files:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        snpfile = open(output_dir + "/snps.txt", 'w')
        swriter = csv.writer(snpfile, delimiter='\t')
        for snp, props in snps.items():
            for frag in props["fragments"]:
                swriter.writerow((snp, props["chr"], props["locus"],
                                  frag["frag"], frag["enzyme"]))
        snpfile.close()

    return snps


def find_interactions(
        snps,
        lib_dir,
        include,
        exclude,
        output_dir,
        suppress_intermediate_files=False):
    """Finds fragments that interact with SNP fragments.

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        hic_data_dir: ../../lib/hic_data
        include: A list of cell lines to query for interactions
        exclude: A list of cell lines to exclude from interaction query
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
        A dict named 'interactions' containing fragments interacting with SNPs
            in each cell line e.g.:
            {'MboI':{
                'rs9462794':{  # SNP rsID
                    'GM12878':{  # Cell line
                        ('6', '30099'): {
                            },
                        ('6', '31188'): [replicate1, replicate2] #Fragment chr and ID
                        }
                    'HeLa':{
                        ('12', '19151'): [replicate1]
                        }
                    }
                'rs12198798':{..}
                'rs6909834':{...}
                }
            'HindIII'{...}
            }

        If suppress_intermediate_files=False, a snp-gene_interaction.txt file is written
          with the ff columns:
            1. SNP rsID
            2. Cell line
            3. Fragment chromosome
            4. Fragment ID
    """
    print("Finding interactions...")
    # Look for all interactions involving SNP fragments in the HiC databases
    if include:
        include_set = set(include)
    if exclude:
        exclude_set = set(exclude)
    interactions = {}
    for enzyme in restriction_enzymes:
        interactions[enzyme] = {}
        for snp in snps:
            interactions[enzyme][snp] = {}
    for enzyme in restriction_enzymes:
        print("\tSearching HiC libraries restricted with " + enzyme)
        hic_data_dir = os.path.join(lib_dir,enzyme, 'hic_data')
        for cell_line in os.listdir(hic_data_dir):
            if (include and cell_line not in include) or (
                    exclude and cell_line in exclude):
                continue  # Enforce the -n or -x options
            if os.path.isdir(os.path.join(hic_data_dir, cell_line)):
                print("\t\tSearching cell line " + cell_line)
                for snp in snps:
                   # A set of unique interactions for SNP
                    for fragment in snps[snp]['fragments']:
                        if not fragment['enzyme'] == enzyme:
                            continue
                        interactions[enzyme][snp][cell_line] = {}
                        print("\t\t\tFinding interactions for " + snp)
                        for replicate in os.listdir(hic_data_dir + '/' + cell_line):
                            if replicate.endswith(".db"):
                                rep_db = sqlite3.connect(
                                    hic_data_dir + '/' + cell_line + '/' + replicate)
                                rep_db.text_factory = str
                                rep_ints = rep_db.cursor()
                                print("\t\t\t\tSearching replicate " + replicate)
                                # Search database for fragments interacting with SNP
                                # fragment

                                from_db = rep_ints.execute(
                                        "SELECT chr2," +
                                        "fragment2 FROM interactions WHERE chr1=?" +
                                        "AND fragment1=?",
                                        [
                                            snps[snp]["chr"],
                                            fragment["frag"]])
                                if from_db:
                                    #interactions[enzyme][snp][cell_line]['replicates'] += 1
                                    for interaction in from_db:
                                        if interaction not in interactions[enzyme][snp][cell_line]:
                                            interactions[enzyme][snp][cell_line]\
                                                [interaction] = set([])
                                        interactions[enzyme][snp][cell_line]\
                                            [interaction].add(replicate)
    if not suppress_intermediate_files:
        intfile = open(output_dir + "/snp-gene_interactions.txt", 'w')
        iwriter = csv.writer(intfile, delimiter='\t')
        for enzyme in interactions:
            for snp in interactions[enzyme]:
                for cell_line in interactions[enzyme][snp]:
                    for interaction in interactions[enzyme][snp][cell_line]:
                        iwriter.writerow(
                            (snp, cell_line, interaction[0], interaction[1],
                             len(interactions[enzyme][snp][cell_line][interaction]), enzyme))
        intfile.close()
    return interactions


def find_genes(
        interactions,
        fragment_database_fp,
        gene_bed_fp,
        output_dir,
        suppress_intermediate_files=False):
    """Identifies genes in fragments that interact with SNP fragments.

    Args:
        interactions: The dictionary fragements that interact with SNP fragments
            returned from find_interactions
        fragment_database_fp: ../../lib/Homo_sapiens.GRCh37.75.dna.fragments.db
        gene_bed_fp: ../../lib/gene_reference.bed
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
        A dict named 'genes' containing genes in fragments that interact with SNP
                fragments e.g.
            {'rs9462794':{  # SNP rsID
                'PHACTR1':{  # Gene
                    'GM12878_Rao2014': {'interactions': 44, 'replicates': 22, 'rep_present': 10},
                    'KBM7_Rao2014': {'interactions': 9, 'replicates': 5, 'rep_present': 2}
                     }
                'EDN1':{...}
                }
             'rs12198798':{..}
             'rs6909834':{...}
             }

        If suppress_intermediate_files=False, SNP-gene pairs with HiC contacts
          in only one libary replicate in only one cell line are written to
          genes_to_remove.txt.The others interactions are written to genes.txt.
          Each file has the ff columns:
            1. SNP rsID
            2. Gene name
            3. Cell line
            4. HiC contact counts
            5. Number of replicates in which contact is found
            6. Number of cell line replicates
    """
    print("Identifying interactions with genes...")
    hs_gene_bed = pybedtools.BedTool(gene_bed_fp)
    genes = {}
    for enzyme in interactions:
        print("\tin libraries restricted with " + enzyme)
        fragment_database_fp = os.path.join(lib_dir, enzyme, 'dna.fragments.db')
        fragment_index_db = sqlite3.connect(fragment_database_fp)
        fragment_index_db.text_factory = str
        fragment_index = fragment_index_db.cursor()
        for snp in interactions[enzyme]:
            print("\t\t"+snp)
                # Generate BED file of all fragments interacting with
                # SNP-containing fragment
            for cell_line in interactions[enzyme][snp]:
                snpgenes_exist = False
                temp_snp_bed = open(os.path.join(output_dir + "/temp_snp_bed.bed"), 'w')
                twriter = csv.writer(temp_snp_bed, delimiter = '\t')
                num_reps = len([rep for rep in os.listdir(
                    os.path.join(lib_dir, enzyme, 'hic_data', cell_line))
                    if rep.endswith('.db')])
                for interaction in interactions[enzyme][snp][cell_line]:
                    fragment_index.execute(
                        "SELECT start, end FROM fragments WHERE chr=? and fragment=?",
                        [interaction[0], interaction[1]])
                    fragment_pos = fragment_index.fetchall()
                    if fragment_pos is None:
                        print(
                            "\tWarning: error retrieving fragment %s on chromosome %s"
                            % (interaction[1], interaction[0]))
                        continue
                    for f in fragment_pos:
                        twriter.writerow(("chr" + interaction[0], f[0], f[1],
                                          interactions[enzyme][snp][cell_line][interaction]))
                    if not snpgenes_exist:
                        snpgenes_exist = True
                temp_snp_bed.close()
                if snpgenes_exist:
                    if not snp in genes:
                        genes[snp] = {}
                    int_bed = pybedtools.BedTool(output_dir + "/temp_snp_bed.bed")
                    # Return a list of genes with which SNP is interacting
                    #    and the number of HiC contacts for each cell line.
                    gene_bed = int_bed.intersect(hs_gene_bed, loj=True)
                    for feat in gene_bed:
                        gene_name = feat[7]
                        if gene_name == '.' or feat[4] == '.' or \
                           feat[5] == '-1' or feat[6] == '-1': # '.' indicates a NULL overlap.
                            continue
                        if not gene_name in genes[snp]:
                            genes[snp][gene_name] = {}
                        if not cell_line in genes[snp][gene_name]:
                            genes[snp][gene_name][cell_line] = {}
                            genes[snp][gene_name][cell_line]['interactions'] = 0
                            genes[snp][gene_name][cell_line]['rep_present'] = []
                        genes[snp][gene_name][cell_line]['interactions'] += 1
                        genes[snp][gene_name][cell_line]['replicates'] = num_reps
                        rep_present = feat[3].replace('}', '')
                        rep_present = rep_present.replace('{', '')
                        rep_present = rep_present.replace("'", "")
                        rep_present = rep_present.split(',')
                        rep_present = [e.strip() for e in rep_present]
                        genes[snp][gene_name][cell_line]['rep_present'] += rep_present
    os.remove(output_dir + "/temp_snp_bed.bed")
    snps_to_remove = {}
    for enzyme in interactions:
        snps_to_remove[enzyme] = []
        for snp in interactions[enzyme]:
            if not snp in genes:
                print("\tNo SNP-gene spatial interactions detected for %s, \
                      removing from analysis" % (snp,))
                snps_to_remove[enzyme].append(snp)
    for enzyme in snps_to_remove: #Update snps and interactions mappings.
        for snp in snps_to_remove[enzyme]:
            for i, frag in enumerate(snps[snp]['fragments']):
                if frag['enzyme'] == enzyme:
                    snps[snp]['fragments'].remove(snps[snp]['fragments'][i])
            del interactions[enzyme][snp]
    genes_to_remove = []
    del_genefile = open(os.path.join(output_dir, 'genes_removed.txt'), 'w')
    dwriter = csv.writer(del_genefile, delimiter = '\t')
    for snp in genes:
        for gene in genes[snp]:
            num_cell_line = len(genes[snp][gene])
            for cell_line in genes[snp][gene]:
                rep_present = len(set(genes[snp][gene][cell_line]['rep_present']))
                interactions = genes[snp][gene][cell_line]['interactions']
                replicates = genes[snp][gene][cell_line]['replicates']
                if interactions/replicates <=1 and rep_present < 2 and\
                num_cell_line < 2:
                    genes_to_remove.append((snp,gene))
                    dwriter.writerow((snp, gene, cell_line, interactions,
                                     rep_present, replicates))
    del_genefile.close()
    for pair in genes_to_remove:
        del genes[pair[0]][pair[1]]
    if not suppress_intermediate_files:
        genefile = open(output_dir + "/genes.txt", 'w')
        gwriter = csv.writer(genefile, delimiter = '\t')
        for snp in genes:
            for gene in genes[snp]:
                for cell_line in genes[snp][gene]:
                    gwriter.writerow((snp, gene, cell_line,
                                      genes[snp][gene][cell_line]['interactions'],
                                      genes[snp][gene][cell_line]['rep_present'],
                                      genes[snp][gene][cell_line]['replicates']))
    return genes

def find_eqtls(
        snps,
        genes,
        eqtl_data_dir,
        gene_database_fp,
        fdr_threshold,
        local_databases_only,
        num_processes,
        output_dir,
        gene_dict_fp,
        snp_dict_fp,
        suppress_intermediate_files=False):
    """Identifies genes in fragments that interact with SNP fragments.

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        genes: The dictionary of genes that interact with SNPs from find_genes
        eqtl_data_dir: ../../eQTLs  #For local eQTL query
        gene_database_fp: ../../gene_reference.db
        fdr_threshold: Significance threshold for the FDR. Default = 0.05
        local_databases_only: boolean specifying if only local eQTL database
            should be querried.
        num_processes: Number of processors to be used when multiprocessing.
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', eqtls.txt file is written to output_dir

    Returns:
        num_sig: Number of eQTL interactions with adj_p_values >= FDR threshold

        If suppress_intermediate_files=False, a eqtls.txt file is written
          with the ff columns:
            1. SNP rsID
            2. Gene name
            3. Tissue of eQTL interaction
            4. eQTL pvalue
            5. eQTL effect size
    """
    print("Identifying eQTLs of interest...")
    eqtls = {}  # A mapping of SNP-gene eQTL relationships
    p_values = []  # A sorted list of all p-values for computing FDR
    num_tests = 0  # Total number of tests done
    try:
        os.remove(os.path.join(output_dir, 'eqtls.txt'))
    except OSError:
        pass
    num_tests, to_online_query = query_local_databases(
        eqtl_data_dir,
        genes,
        gene_dict_fp,
        snp_dict_fp,
        p_values,
        output_dir)
    print(num_tests, len(p_values), len(to_online_query))
    if not local_databases_only:
        num_tests += query_GTEx_service(snps,
                                        genes,
                                        to_online_query,
                                        p_values,
                                        num_processes,
                                        output_dir)
    
    return num_tests, p_values

def query_local_databases(
        eqtl_data_dir, genes, gene_dict_fp, snp_dict_fp, p_values, output_dir):
    """Not used at the moment.
    """
    print("\tQuerying local databases.")
    num_tests = 0
    to_online_query = []
    gene_dict_db = sqlite3.connect(gene_dict_fp)
    gene_dict_db.text_factory = str
    gene_dict = gene_dict_db.cursor()
    snp_dict_db = sqlite3.connect(snp_dict_fp)
    snp_dict_db.text_factory = str
    snp_dict = snp_dict_db.cursor()
    eqtlfile = open(os.path.join(output_dir, 'eqtls.txt'), 'a')
    ewriter = csv.writer(eqtlfile, delimiter='\t')
    for db in os.listdir(
            eqtl_data_dir):  # Iterate through databases of eQTLs by tissue type
        tissue = db[:db.rfind('.')]
        print("\t\tQuerying " + tissue)
        eqtl_index_db = sqlite3.connect(os.path.join(eqtl_data_dir, db))
        eqtl_index_db.text_factory = str
        eqtl_index = eqtl_index_db.cursor()
        for snp in genes:
            for gene in genes[snp]:
                try:
                    snp_dict.execute("SELECT variant_id FROM lookup WHERE rsID = ?;",
                                     (snp,))
                    variant_id = snp_dict.fetchone()[0]
                    gene_dict.execute("SELECT gene_id FROM genes WHERE gene_name = ?;",
                                      (gene,))
                    gene_id = gene_dict.fetchone()[0]
                    eqtl_index.execute(
                        "SELECT pval_nominal, slope FROM associations " +\
                        "WHERE variant_id=? AND gene_id=?",
                        (variant_id, gene_id))
                    if not eqtl_index.fetchall():
                        to_online_query.append((snp, gene, tissue))
                        continue
                    for eqtl in eqtl_index.execute(
                            "SELECT pval_nominal, slope FROM associations " +\
                            "WHERE variant_id=? AND gene_id=?",
                            (variant_id, gene_id)):
                        ewriter.writerow((
                            snp, gene, tissue, eqtl[0], eqtl[1]))
                        bisect.insort(p_values, eqtl[0])
                        num_tests += 1
                except:
                    to_online_query.append((snp, gene, tissue))
        eqtl_index.close()
        eqtl_index_db.close()
    eqtlfile.close()
    gene_dict.close()
    gene_dict_db.close()
    snp_dict.close()
    snp_dict_db.close()
    return num_tests, to_online_query

def query_GTEx_service(
        snps,
        genes,
        to_query,
        p_values,
        num_processes,
        output_dir):
    """Queries GTEx for eQTL association between SNP and gene.

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        genes: The dictionary of genes that interact with SNPs from find_genes
        eqtls: A mapping of eQTl relationships between SNPs to genes, which
            further  maps to a list of tissues in which these eQTLs occur
        p_values: A sorted list of all p-values for use computing FDR.
        num_processes: Number of processors to be used when multiprocessing.
        output_dir: User-specified directory for results. Defaults to inputs
            directory.

    Returns:
        num_tests: The number of successful GTEx responses
        eqtls and p_values are updated.
        If GTEx query fails for SNP-gene-tissue request, the failed requests
            are written to failed_GTEx_requests.txt thus:
            [
                {
                    'variantId': 'rs9462794',
                    'gencodeId': LINC00669'
                    'tissueSiteDetailId': 'Artery_Aorta'
                }
            ]
    """

    if num_processes > 2 : # Ensure GTEx is not 'flooded' with requests
        num_processes = 2
    manager = multiprocessing.Manager()
    num_tests = 0
    reqLists = [[]]
    for assoc in to_query:
        if len(reqLists[-1]) < 950: # >1000 requests is buggy.
            reqLists[-1].append({"variantId": assoc[0],
                                 "gencodeId": assoc[1],
                                 "tissueSiteDetailId": assoc[2]})
        elif reqLists[-1][-1]["variantId"] == assoc[0] and \
             reqLists[-1][-1]["gencodeId"] == assoc[1]:
            # Ensure that all tissues of a SNP-gene pair are in same ReqList
            reqLists[-1].append({"variantId": assoc[0],
                                 "gencodeId": assoc[1],
                                 "tissueSiteDetailId": assoc[2]})
        else:
            reqLists.append([{"variantId": assoc[0],
                              "gencodeId": assoc[1],
                              "tissueSiteDetailId": assoc[2]}])
    print("\tQuerying GTEx online service...")
    gtexResponses = manager.list()
    procPool = multiprocessing.Pool(processes=num_processes)
    for i, reqList in enumerate(reqLists, start=1):
        procPool.apply_async(
            send_GTEx_query, (i, len(reqLists), reqList, gtexResponses))
        time.sleep(10)
    procPool.close()
    procPool.join()
    print("\t\tGTEx responses received: " + str(len(gtexResponses)))
    print('Writing eQTL file...')
    eqtlfile = open(output_dir + '/eqtls.txt', 'a')
    ewriter = csv.writer(eqtlfile, delimiter = '\t')
    failed_requests = []
    for response in gtexResponses:
        try:
            results = response[1].json()["result"]
            #results += response[1].json()["result"]
            for result in results:
                if (str(result["geneSymbol"]) == "gene not found" or
                    not result["snpId"] in genes or
                    result["pValue"] == "N/A"):
                    continue
                num_tests += 1
                snp = result["snpId"]
                gene = result["geneSymbol"]
                p = float(result["pValue"])
                effect_size = float(result["nes"])
                bisect.insort(p_values, p)
                ewriter.writerow((snp,
                                  gene,
                                  result["tissueSiteDetailId"],
                                  p,
                                  effect_size))
        except Exception as e:
            print("\t\tWARNING: bad response (%s)" % response[1])
            failed_requests.append(response[0])
    eqtlfile.close()
    print('\n')    
    if failed_requests:
        # TODO: Handle the failed requests procedure
        with open(output_dir + "/failed_GTEx_requests.txt", 'w') \
             as failed_requests_file:
            failed_requests_file.write(str(failed_requests) + '\n')
    return num_tests

def calc_hic_contacts(snp_gene_dict):
    """Calculates score of HiC contacts between SNP and gene.

    Args:
        genes: Dictionary from find_genes
        snp:
        gene:

    Returns:
        hic_score: sum of the averages of contacts per cell line.
        score_list: a list of the means of contacts in each cell line
    """
    hic_score = 0
    cell_lines = sorted(snp_gene_dict.keys())
    scores = []
    for cell_line in cell_lines:
        score = snp_gene_dict[cell_line]['interactions'] / float(snp_gene_dict[
                                                    cell_line]['replicates'])
        scores.append('{:.2f}'.format(score))
        hic_score += score
    return (', '.join(cell_lines), ', '.join(scores), 
            "{:.4f}".format(hic_score))

def send_GTEx_query(num, num_reqLists, reqList, gtexResponses):
    """Posts and receives requests from GTEx

    Args:
        num: An integer indicating the request being processed.
        num_reqLists: The total number of requests to be submitted
        reqList: A list of < 950 SNP, gene and tissue JSON data
        gtexResponses: A manager.list() object to handle requests.

    Returns:
        gtexResponses:
           [
              (ReqList:[]
               {'beta': '0.0832273122097661',
                'gencodeId': 'ENSG00000137872.11',
                'geneSymbol': 'SEMA6D',
                'pvalue': '0.37233196039320215',
                'pvalueThreshold': None,
                'se': '0.0929115299234',
                'variantId': 'rs9462794',
                'tissueId': 'Prostate',:
                'tstat': '0.8957694731579682',
                'variantId': '6_12445847_A_T_b37'})
           ]
    """
    s = requests.Session()
    s.verify = GTEX_CERT
    try:
        while True:
            print("\t\tSending request %s of %s" % (num, num_reqLists))
            gtex_url = "https://gtexportal.org/api/b0.9/association/dyneqtl"
            res = s.post(gtex_url, json=reqList)
            #print("Printing status code...", res.status_code )
            if res.status_code == 200:
                gtexResponses.append((reqList, res))
                time.sleep(60)
                return
            elif res.status_code == 500 or res.status_code == 400:
                print("\t\tThere was an error processing request %s. \
                      Writing to failed request log and continuing." % num)
                gtexResponses.append((reqList, "Processing error"))
                time.sleep(60)
                return
            else:
                print("\t\tRequest %s received response with status %s. "+ \
                      "Trying again in five minutes." % (num, res.status_code))
                time.sleep(300)
    except requests.exceptions.ConnectionError:
        try:
            print("\t\tWarning: Request %s experienced a connection error. \
                  Retrying in five minutes." % num)
            time.sleep(300)
            while True:
                print("\t\tSending request %s of %s" % (num, num_reqLists))
                res = s.post(
                    "http://gtexportal.org/api/v6p/dyneqtl?v=clversion",
                    json=reqList)
                if res.status_code == 200:
                    gtexResponses.append((reqList, res))
                    time.sleep(30)
                    return
                elif res.status_code == 500:
                    print("\t\tThere was an error processing request %s. \
                          Writing to failed request log and continuing." % num)
                    gtexResponses.append((reqList, "Processing error"))
                    time.sleep(60)
                    return
                else:
                    print("\t\tRequest %s received response with status: %s. \
                          Writing to failes request log and continuing." \
                          % (num, res.status_code))
                    gtexResponses.append((reqList, res.status_code))
                    time.sleep(60)
                    return
        except requests.exceptions.ConnectionError:
            print("\t\tRetry failed. Continuing, but results will be incomplete.")
            gtexResponses.append((reqList, "Connection failure"))
            time.sleep(60)
            return
    s.close()

def get_gene_expression_information(eqtls, expression_table_fp, output_dir):

    print("Getting gene expression information...")
    gene_df = pandas.read_table(expression_table_fp, index_col='Symbol')
    gene_exp = pandas.DataFrame(data=None, columns=gene_df.columns)
    for snp in eqtls:
        for gene in eqtls[snp]:
            if gene not in gene_exp.index:
                try:
                    gene_exp = gene_exp.append(gene_df.ix[gene])
                except KeyError:
                    print("Warning: No gene expression information for %s" 
                          % gene)
    gene_exp.to_csv(
        path_or_buf=output_dir +
        "/gene_expression_table.txt",
        sep='\t')

def get_gene_expression_extremes(gene, gene_exp):
    try:
        max_expression = max(gene_exp.iteritems(), 
                             key=lambda exp: max(exp[1]))
    except KeyError:
        return tuple(['NA', 'NA', 'NA', 'NA'])
    min_expression = min(gene_exp.iteritems(), 
                         key=lambda exp: min(exp[1]))
    return tuple([max_expression[0], max(max_expression[1]), min_expression[0], 
                 min(min_expression[1])])

def get_tissue_expression(gene_tissue_tpl):
    try:
        return list(GENE_DF.at[gene_tissue_tpl[0], gene_tissue_tpl[1]])
    except TypeError: 
        return  list([GENE_DF.at[gene_tissue_tpl[0], gene_tissue_tpl[1]]])
    except KeyError:
        print("\t\tWarning: No expression information for %s in %s." % (gene, 
                                                                       tissue))
        return None 
    
def produce_summary(
        p_values, snps, genes, gene_database_fp,
        expression_table_fp, fdr_threshold, output_dir, buffer_size_in,
        buffer_size_out, num_processes):    
    """Write final results of eQTL-eGene associations

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        genes: The dictionary of genes that interact with SNPs from find_genes
        gene_database_fp: ../../gene_reference.db
        expression_table_fp: ../../GTEx
        output_dir:

    Returns:
        summary.txt: A file with the ff structure.
            1. SNP
            2. SNP chromosome
            3. SNP locus
            4. Gene name
            5. Gene chromosome
            6. Gene start
            7. Gene end
            8. Tissue
            9. p-value
            10. Adj p-value
            11. Effect size
            12. Cell lines
            13. GTEx_cis_p_Threshold
            14. cis_SNP-gene_interaction
            15. SNP-gene_Distance
            16. Expression_Level_In_eQTL_Tissue
            17. Max_Expressed_Tissue
            18. Maximum_Expression_Level
            19. Min_Expressed_Tissue
            20. Min_Expression_Level

        sig_eqtls.text: A file with eQTL associations with
            adj_p_values <= FDR threshold. Same structure as summary.txt
    """
    num_sig = {} 
    # Number of eQTLs deemed significant under the given threshold
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    summary = open(output_dir + "/summary.txt", 'w', buffering=buffer_size_out)
    summ_writer = csv.writer(summary, delimiter='\t')
    summ_header = ['SNP',
                   'SNP_Chromosome',
                   'SNP_Locus',
                   'Gene_Name',
                   'Gene_Chromosome',
                   'Gene_Start',
                   'Gene_End',
                   'Tissue',
                   'p-value',
                   'Adj_p-value',
                   'Effect_Size',
                   'Cell_Lines',
                   'Cell_Lines_HiC_scores',
                   'HiC_Contact_Score',
                   'GTEx_cis_p_Threshold',
                   'cis_SNP-gene_interaction',
                   'SNP-gene_Distance',
                   'Expression_Level_In_eQTL_Tissue',
                   'Max_Expressed_Tissue',
                   'Maximum_Expression_Level',
                   'Min_Expressed_Tissue',
                   'Min_Expression_Level']
    summ_writer.writerow(summ_header)
    sig_file = open(os.path.join(output_dir, 'significant_eqtls.txt'), 'w', 
                    buffering=buffer_size_out)
    sigwriter = csv.writer(sig_file, delimiter='\t')
    sigwriter.writerow(summ_header) 
    gene_index_db = sqlite3.connect(gene_database_fp)
    gene_index_db.text_factory = str
    col_names = [res[1]
                 for res in gene_index_db.execute("PRAGMA table_info(genes)")]
    include_p_thresh = "p_thresh" in col_names
    if include_p_thresh:
        query_str = ("SELECT chr, start, end, p_thresh FROM genes " 
                     "WHERE symbol=?")
    else:
        query_str = "SELECT chr, start, end FROM genes WHERE symbol=?"
    gene_exp = {} 
    eqtlfile = open(os.path.join(output_dir, 'eqtls.txt'), 'r', 
                    buffering=buffer_size_in)
    ereader = csv.reader(eqtlfile, delimiter = '\t')
    
    p_values_map = {} 
    print('Adjusting p-values...')
    adjusted_p_values = compute_adj_pvalues(p_values)
    for i in range(len(p_values)):
        p_values_map[p_values[i]] = adjusted_p_values[i]

    to_file = [] 
    genes_from_file = []
    snps_from_file = [] 
    
    for i, row in enumerate(ereader):
        snp = row[0]
        gene = row[1]
        tissue = row[2]
        p_val = row[3]
        effect_size = row[4]
        genes_from_file.append(gene)
        snps_from_file.append(snp)
        qvalue = p_values_map[float(p_val)]
        snp_info = {
            "chr": snps[snp]["chr"],
            "locus": snps[snp]["locus"],
            "fragments": snps[snp]["fragments"]}
        gene_chr = "NA"
        gene_start = "NA"
        gene_end = "NA"
        max_length = 0
 
        gene_stat = gene_index_db.execute(query_str, (gene,)).fetchone()
        # Consider "canonical" to be the longest record where
        # multiple records are present
        if abs(gene_stat[2] - gene_stat[1]) > max_length:
            if gene_stat[0].startswith("chr"):
                gene_chr = gene_stat[0][3:]
            else:
                gene_chr = gene_stat[0]
            gene_start = gene_stat[1]
            gene_end = gene_stat[2]
        max_length = gene_stat[2] - gene_stat[1]
        if include_p_thresh:
            p_thresh = gene_stat[3]
        else:
            p_thresh = "NA"
        # eQTL is cis if the SNP is within 1Mbp of the gene
        cis = gene_chr == snp_info["chr"] and (
            (snp_info["locus"] > gene_start - 1000000) and
            (snp_info["locus"] < gene_end + 1000000))
        cell_lines = ''
        try:
            cell_lines = list(genes[snp][gene])
        except KeyError:
            cell_lines = "NA"
        p_thresh = p_thresh
        interaction = True
        if cis:
            interaction = True
        else:
            interaction = False
        distance_from_snp = 0
        if(gene_chr == "NA"):
            distance_from_snp = "NA"  # If gene location information is missing
        elif(not snp_info["chr"] == gene_chr):
            distance_from_snp = "NA"  # Not applicable to trans interactions
        elif(snp_info["locus"] < gene_start):
            distance_from_snp = gene_start - snp_info["locus"]
        elif(snp_info["locus"] > gene_end):
            distance_from_snp = snp_info["locus"] - gene_end
        if not gene in gene_exp:
            gene_exp[gene] = {} 
        if not tissue in gene_exp[gene]:
            gene_exp[gene][tissue] = 'NA'
        to_file.append([snp,
                   snp_info["chr"],
                   snp_info["locus"],
                   gene,
                   gene_chr,
                   gene_start,
                   gene_end,
                   tissue,
                   p_val,
                   qvalue,
                   effect_size,
                   p_thresh,
                   interaction,
                   distance_from_snp])

    
    genes_from_file = dict.fromkeys(genes_from_file).keys()
    snps_from_file = dict.fromkeys(snps_from_file).keys()
    snps_genes = [(snp, gene) for snp in snps_from_file 
                    for gene in genes_from_file]
    global GENE_DF
    # Efficiently allowing access to a shared namespace by individual processes
    GENE_DF = pandas.read_table(expression_table_fp, index_col="Description", 
                                 engine='c', compression=None, memory_map=True)
    all_tissues = list(GENE_DF) 
    genes_tissues = [(gene, tissue) for gene in genes_from_file 
                        for tissue in all_tissues]
    
    current_proc = psutil.Process()
    pool = multiprocessing.Pool(processes=min(num_processes,
                                            len(current_proc.cpu_affinity())))

    print("Collecting gene expression rates...")
    gene_tissue_expression_extremes = pool.map(get_tissue_expression, 
                                             genes_tissues)
    for i in range(len(genes_tissues)):
        expression_rates = [x for x in gene_tissue_expression_extremes[i] 
                            if x > 0.0]
        if expression_rates:
            gene_exp[genes_tissues[i][0]][
                     genes_tissues[i][1]] = expression_rates 
    del GENE_DF
    
    print("Determining gene expression extremes...")
    for gene in genes_from_file:
        gene_exp[gene]['extremes'] = get_gene_expression_extremes(gene,
                                                                gene_exp[gene])  

    print("Computing HiC data...")
    hic_data = pool.map(calc_hic_contacts, 
                      [genes[snp][gene] for snp, gene in snps_genes])
    hic_dict = {} 
    for i in range(len(snps_genes)):
        hic_dict[snps_genes[i]] = hic_data[i]

    print("Writing to summary files...")
    for line in to_file:
        snp = line[0]
        gene = line[3]
        tissue = line[7]
        line.insert(11, hic_dict[(snp, gene)][2])
        line.insert(11, hic_dict[(snp, gene)][1])
        line.insert(11, hic_dict[(snp, gene)][0])
        try:
            line.append(max(gene_exp[gene][tissue]))
        except KeyError:
            line.append('NA')
        line.extend(gene_exp[gene]['extremes'])
        summ_writer.writerow(line)
        if line[9] < fdr_threshold:
            try:
                num_sig[snp] += 1
            except KeyError:
                num_sig[snp] = 1
            sigwriter.writerow(line)
    summary.close()
    sig_file.close()

    return (num_sig)

    
def compute_adj_pvalues(p_values):
    """ A Benjamini-Hochberg adjustment of p values of SNP-gene eQTL
           interactions from GTEx.

    Args:
        p_values: List of sorted p values of all eQTL associations

    Returns:
        adj_pvalues: A corresponding list of adjusted p values to p_values.
    """
    return R.r['p.adjust'](p_values, method='BH')


def produce_overview(genes, eqtls, num_sig, output_dir):
    #TODO: Make compatible without 'eqtls'
    """Generates overview graphs and table

    """
    print("\nProducing overview...")
    stat_table = open(output_dir + "/overview.txt", 'w')
    stat_table.write(
        "SNP\tChromosome\tLocus\tTotal_SNP-gene_Pairs\tTotal_eQTLs\n")
    for snp in eqtls:
        stat_table.write(snp +
                         '\t' +
                         eqtls[snp]["snp_info"]["chr"] +
                         '\t' +
                         str(eqtls[snp]["snp_info"]["locus"]) +
                         '\t' +
                         str(len(genes[snp])) +
                         '\t' +
                         str(len(eqtls[snp]) -
                             1) +
                         '\n')
    stat_table.close()
    print("\tProducing graphs")
    style.use("ggplot")
    if not os.path.isdir(output_dir + "/plots"):
        os.mkdir(output_dir + "/plots")
    int_colours = "rgb"
    eqtl_colours = "myc"
    int_colours = cycle(int_colours)
    eqtl_colours = cycle(eqtl_colours)
    snps_by_chr = {}
    for snp in eqtls:
        chrm = eqtls[snp]["snp_info"]["chr"]
        if not chrm in snps_by_chr:
            snps_by_chr[chrm] = []
        snps_by_chr[chrm].append((eqtls[snp]["snp_info"]["locus"], snp))

    int_colour_list = []
    eqtl_colour_list = []
    num_snpgenes = []
    num_eqtls = []
    rsIDs = []

    chrs = []
    for c in snps_by_chr:
        chrs.append(c)
    # So that the chromosomes are in a logical order on the graph
    chrs.sort(key=natural_keys)
    chr_locs = []
    chr_ticks = []
    chrm_pos = 0
    last_count = 0
    count = 0

    for chrm in chrs:
        snp_list = snps_by_chr[chrm]
        snp_list.sort()  # Sort by locus
        int_colour = next(int_colours)
        eqtl_colour = next(eqtl_colours)
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
    plt.bar(range(len(num_snpgenes)), num_snpgenes, color=int_colour_list)
    plt.bar(range(len(num_eqtls)), num_eqtls, color=eqtl_colour_list)
    axes = plt.gca()
    # So that eQTL values won't appear negative on plot
    axes.yaxis.set_major_formatter(FuncFormatter(abs_value_ticks))
    plt.vlines(chr_locs, axes.get_ylim()[0], axes.get_ylim()[1], colors="k")
    axes.set_xticks(range(len(rsIDs)))
    axes.set_xticklabels(rsIDs, rotation='vertical')
    ax2 = axes.twiny()
    ax2.set_xlim(axes.get_xlim())
    ax2.set_xticks(chr_locs)
    ax2.set_xticklabels(chr_ticks)
    plt.tight_layout()
    plt.savefig(
        output_dir +
        "/plots/snpgene_and_eqtl_overview.png",
        dpi=300,
        format="png")
    plt.clf()


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split('(\d+)', text)]


def abs_value_ticks(x, pos):
    return abs(x)


def retrieve_pathways(eqtls, fdr_threshold, num_processes, output_dir):
    print("Retrieving pathway information...")
    manager = multiprocessing.Manager()
    procPool = multiprocessing.Pool(processes=num_processes)
    pathways = {}  # Keeps track of all pathways covered by genes with a statistically signficant eQTL relationship with a query SNP
    pwresults = manager.list()

    for snp in eqtls:
        for gene in eqtls[snp]:
            if not gene == "snp_info":
                for tissue in eqtls[snp][gene]["tissues"]:
                    if eqtls[snp][gene]["tissues"][tissue]["qvalue"] < fdr_threshold:
                        procPool.apply_async(
                            get_wikipathways_response, [
                                snp, gene, tissue, pwresults])
    procPool.close()
    procPool.join()

    for response in pwresults:
        snp = response[0]
        gene = response[1]
        tissue = response[2]
        for pwres in response[3]:
            pwid = pwres["identifier"]
            if not pwid in pathways:
                pathways[pwid] = {}
                pathways[pwid]["name"] = pwres["name"]
                pathways[pwid]["genes"] = {}
            if not gene in pathways[pwid]["genes"]:
                pathways[pwid]["genes"][gene] = {}
            if not snp in pathways[pwid]["genes"][gene]:
                pathways[pwid]["genes"][gene][snp] = set([])
            pathways[pwid]["genes"][gene][snp].add(tissue)

    with open(output_dir + "/pathways.txt", 'w') as pwfile:
        pwfile.write(
            "WikiPathways_ID\tPathway_Name\tGene_Symbol\tSNP\tTissue\n")
        for pwid in pathways:
            for gene in pathways[pwid]["genes"]:
                for snp in pathways[pwid]["genes"][gene]:
                    for tissue in pathways[pwid]["genes"][gene][snp]:
                        pwfile.write(
                            "%s\t%s\t%s\t%s\t%s\n" %
                            (pwid, pathways[pwid]["name"], gene, snp, tissue))
    return pathways


def get_wikipathways_response(snp, gene, tissue, pwresults):
    wikipathways = wikipathways_api_client.WikipathwaysApiClient()
    kwargs = {
        'query': gene,
        'organism': 'http://identifiers.org/taxonomy/9606'  # Homo sapiens
    }
    res = wikipathways.find_pathways_by_text(**kwargs)
    if res:
        pwresults.append([snp, gene, tissue, res])


def parse_snps_files(snps_files):
    snps = {}
    for snp_file in snps_files:
        with open(snp_file, 'r') as snpfile:
            for line in snpfile:
                snp = line.strip().split('\t')
                snps[snp[0]] = {"chr": snp[1],
                                "locus": int(snp[2]), "frag": snp[3]}
    return snps


def parse_interactions_files(interactions_files):
    interactions = {}
    for interactions_file in interactions_files:
        with open(interactions_file, 'r') as intfile:
            for line in intfile:
                interaction = line.strip().split('\t')
                snp = interaction[0]
                cell_line = interaction[1]
                if not snp in interactions:
                    interactions[snp] = {}
                if not cell_line in interactions[snp]:
                    interactions[snp][cell_line] = set([])
                interactions[snp][cell_line].add(
                    (interaction[2], int(interaction[3])))
    return interactions


def parse_genes_files(genes_files):
    genes = {}
    for gene_file in genes_files.split(' '):
        genefile = open(gene_file, 'r')
        greader = csv.reader(genefile, delimiter = '\t')
        for row in greader:
            gene = row[1]
            snp = row[0]
            cell_line = row[2]
            interactions = int(row[3])
            interactions_replicates = len(row[4].split())
            replicates_count = int(row[5])
            try:
                genes[snp][gene][cell_line] = {
                    'interactions': interactions,
                    'replicates': replicates_count,
                    'rep_present': interactions_replicates} 
            except KeyError:
                if not snp in genes:
                    genes[snp] = {}
                if not gene in genes[snp]:
                    genes[snp][gene] = {}
                genes[snp][gene][cell_line] = {
                    'interactions': interactions,
                    'replicates': replicates_count,
                    'rep_present': interactions_replicates}

    return genes


def parse_eqtls_files(
        eqtls_files, snp_database_fp, gene_database_fp,
        restriction_enzymes, lib_fp, output_dir, fdr_threshold=0.05):
    print("Parsing eQTLs...")
    eqtls = {}
    snps = {}
    genes = {}
    p_values = []  # A sorted list of all p-values for use computing FDR
    num_tests = 0  # Total number of tests done
    num_sig = {}  # Number of eQTLs deemed significant under the given threshold
    # eqtls_files = eqtls_files.split(' ') #@Tayaza: Separate individual eqtl
    # files
    snp_dp = sqlite3.connect(snp_database_fp)
    snp_dp.text_factory = str
    snp_index = snp_dp.cursor()
    gene_dp = sqlite3.connect(gene_database_fp)
    gene_dp.text_factory = str
    gene_index = snp_dp.cursor()
    file_count = 1
    genes_files = []
    for eqtls_file in eqtls_files:
        file_path = eqtls_file[:len(eqtls_file)-9]        
        print('\tParsing eQTL file {} of {}'.format(file_count, len(eqtls_files)))
        file_count += 1
        if os.path.isfile(os.path.join(file_path, 'genes.txt')):
            #TODO: merge genes from multiple genes.txt files
            genes = parse_genes_files(os.path.join(file_path, 'genes.txt'))
            genes_files.append(os.path.isfile(os.path.join(file_path, 'genes.txt')))
        else:
            print('\t\tOops!: We can\'t find  \'genes.txt\' file in the same' +\
                  ' directory as \'eqtls.txt\'. \n\t\tSome info will be missing.')
            #print('Oops!: Please ensure \'genes.txt\' file is in the same' +\
            #      ' directory as \'eqtls.txt\'\n Program terminating.')
        

        eqtlfile =  open(eqtls_file, 'r')
        ereader = csv.reader(eqtlfile, delimiter='\t')
        # Do line count for progress meter
        for i, row in enumerate(ereader, 1):
            snp = row[0]
            gene = row[1]
            tissue = row[2]
            effect_size = float(row[4])
            bisect.insort(p_values, float(row[3]))
            if not snp in snps:
                if os.path.isfile(os.path.join(file_path, 'snps.txt')):
                    snp_file = open(os.path.join(file_path, 'snps.txt'), 'r')
                    sreader = csv.reader(snp_file, delimiter = '\t')
                    for row in sreader:
                        if row[0] == snp:
                            try:
                                snps[snp]['fragments'].append({
                                    'frag': row[3], 'enzyme': row[4]})
                            except KeyError:
                                snps[snp] = {}
                                snps[snp]['chr'] = row[1]
                                snps[snp]['locus'] = int(row[2])
                                snps[snp]['fragments'] = []
                                snps[snp]['fragments'].append({
                                    'frag': row[3], 'enzyme': row[4]})

                else:  # In case there is no snps.txt file
                    snp_index.execute("SELECT * FROM snps WHERE rsID=?", (snp,))
                    from_db = snp_index.fetchone()
                    snps[snp] = {'chr': from_db[1],
                             'locus': from_db[2],
                             'fragments': []} #Handle fragments from multiple REs
                    enzymes = restriction_enzymes.split(',')
                    enzymes = [e.strip() for e in enzymes]
                    for enzyme in enzymes:
                        conn = sqlite3.connect(os.path.join(
                            lib_fp, enzyme, 'dna.fragments.db'))
                        conn.text_factory = str
                        cur = conn.cursor()
                        cur.execute(
                        "SELECT fragment FROM fragments WHERE chr=? "+\
                            "AND start<=? AND end>=?", (snps[snp]['chr'],
                                                        snps[snp]['locus'],
                                                        snps[snp]['locus']))
                        snps[snp]['fragments'].append({
                            'frag': cur.fetchone()[0],
                            'enzyme': enzyme})
                        cur.close()
                        conn.close()
        eqtlfile.close()
    snp_index.close()
    snp_dp.close()
    #gene_index.close()
    #gene_dp.close()
    joined_eqtl_file = open(os.path.join(output_dir, 'eqtls.txt'), 'w')
    for eqtls_file in eqtls_files:
        with open(eqtls_file, 'r') as efile:
            shutil.copyfileobj(efile, joined_eqtl_file)
    
    return (p_values, snps, genes)


def build_snp_index(
        snp_dir,
        output_fp,
        config,
        id_col=4,
        chr_col=1,
        locus_col=2,
        do_not_tidy_up=False):
    if not output_fp:
        output_fp = config["SNP_DATABASE_FP"]

    print("Building SNP index...")
    if not os.path.isdir(snp_dir):
        print("Error: argument to build SNP index must be a directory.")
        return

    if os.path.isfile(output_fp):
        upsert = input(
            "WARNING: Upserting input to existing SNP database (%s). Continue? [y/N] " %
            output_fp)
        if not upsert.lower() == 'y':
            print("Did not write to existing SNP database.")
            return

    if not os.path.isdir(os.path.dirname(output_fp)):
        os.makedirs(os.path.dirname(output_fp))

    snp_index_db = sqlite3.connect(output_fp)
    snp_index = snp_index_db.cursor()
    snp_index.execute(
        "CREATE TABLE IF NOT EXISTS snps (rsID text unique, chr text, locus integer)")
    snp_index.execute("CREATE INDEX IF NOT EXISTS id ON snps (rsID,chr,locus)")
    res = ""
    id_col = int(id_col) - 1
    chr_col = int(chr_col) - 1
    locus_col = int(locus_col) - 1
    for bed_fp in os.listdir(snp_dir):
        print("\tProcessing " + bed_fp)
        bed = open(snp_dir + '/' + bed_fp, 'r')
        for line in bed:
            if not line.startswith("chr"):
                continue
            snp = line.strip().split('\t')
            try:
                snp_index.execute("INSERT INTO snps VALUES(?,?,?)", [
                                  snp[id_col], snp[chr_col][snp[chr_col].find("chr") + 3:], int(snp[locus_col])])
            except sqlite3.IntegrityError:
                if res == "":
                    res = input(
                        "Warning: database already contains an entry for SNP with ID %s. Overwrite?\n1: Yes\n2: No (default)\n3: Yes To All\n4: No To All\nChoice: " %
                        snp[id_col])
                if res.strip() == "1":
                    print("Overwriting SNP %s" % snp[id_col])
                    snp_index.execute(
                        "DELETE FROM snps WHERE rsID=?", [
                            snp[id_col]])
                    snp_index.execute("INSERT INTO snps VALUES(?,?,?)", [
                                      snp[id_col], snp[chr_col][snp[chr_col].find("chr") + 3:], int(snp[locus_col])])
                    res = ""
                elif res.strip() == "3":
                    print("Overwriting SNP %s" % snp[id_col])
                    snp_index.execute(
                        "DELETE FROM snps WHERE rsID=?", [
                            snp[id_col]])
                    snp_index.execute("INSERT INTO snps VALUES(?,?,?)", [
                                      snp[id_col], snp[chr_col][snp[chr_col].find("chr") + 3:], int(snp[locus_col])])
                elif res.strip() == "4":
                    print("Skipping input SNP %s" % snp[id_col])
                    pass
                else:
                    print("Skipping input SNP %s" % snp[id_col])
                    res = ""
        bed.close()
    print("\tWriting SNP index to file...")
    snp_index_db.commit()
    print("Done building SNP index.")
    if not do_not_tidy_up:
        print("Tidying up...")
        shutil.rmtree(snp_dir)


def build_hic_index(
        input_hic_fp,
        output_fp=None,
        chr1_col=3,
        chr2_col=7,
        frag1_col=5,
        frag2_col=9,
        mapq1_col=10,
        mapq2_col=11,
        mapq_cutoff=30,
        do_not_tidy_up=False):
    if not output_fp:
        if not input_hic_fp.rfind('.') == -1:
            output_fp = input_hic_fp[:input_hic_fp.rfind('.')] + ".db"
        else:
            output_fp = input_hic_fp + ".db"

    if os.path.isfile(output_fp):
        overwrite = input(
            "WARNING: Overwriting existing HiC database (%s). Continue? [y/N] " %
            output_fp)
        if not upsert.lower() == 'y':
            print("Did not overwrite existing HiC database.")
            return
        os.remove(output_fp)

    # Do line count for progress meter
    print("Determining table size...")
    hic_table = open(input_hic_fp, 'r')
    lines = 0
    for i in hic_table:
        lines += 1
    hic_table.close()
    lines = lines // 100 * 100  # Get an approximation
    do_linecount = not lines == 0

    int_db = sqlite3.connect(output_fp)
    interactions = int_db.cursor()
    interactions.execute(
        "CREATE TABLE IF NOT EXISTS interactions (chr1 text, fragment1 text, chr2 text, fragment2 text)")
    interactions.execute(
        "CREATE INDEX IF NOT EXISTS i_1 ON interactions (chr1, fragment1)")
    chr1_col = int(chr1_col) - 1
    chr2_col = int(chr2_col) - 1
    frag1_col = int(frag1_col) - 1
    frag2_col = int(frag2_col) - 1
    mapq1_col = int(mapq1_col) - 1
    mapq2_col = int(mapq2_col) - 1
    with open(input_hic_fp, 'r') as rao_table:
        print("Indexing HiC interaction table...")
        for i, line in enumerate(rao_table):
            if do_linecount:
                if i % (lines / 100) == 0:
                    print("\tProcessed %d%%..." % ((float(i) / float(lines)) * 100))
            interaction = line.strip().split(' ')
            if int(
                    interaction[mapq1_col]) >= mapq_cutoff and int(
                    interaction[mapq2_col]) >= mapq_cutoff:
                chr1 = interaction[chr1_col]
                frag1 = interaction[frag1_col]
                chr2 = interaction[chr2_col]
                frag2 = interaction[frag2_col]
                interactions.execute(
                    "INSERT INTO interactions VALUES(?,?,?,?)", [
                        chr1, frag1, chr2, frag2])
                interactions.execute(
                    "INSERT INTO interactions VALUES(?,?,?,?)", [
                        chr2, frag2, chr1, frag1])
    int_db.commit()
    interactions.close()
    print("Done indexing HiC interaction table.")
    if not do_not_tidy_up:
        print("Tidying up...")
        os.remove(input_hic_fp)


def digest_genome(
        genome,
        restriction_enzyme,
        output_fp,
        output_db,
        do_not_index=False,
        linear=False):
    if not output_fp:
        if not genome.rfind('.') == -1:
            output_fp = "%s.%s.fragments.bed" % (
                genome[:genome.rfind('.')], restriction_enzyme)
        else:
            output_fp = "%s.%s.fragments.bed" % (genome, restriction_enzyme)
    if not os.path.isdir(os.path.dirname(output_fp)):
        os.makedirs(os.path.dirname(output_fp))
    if os.path.isfile(output_fp):
        overwrite = input(
            "WARNING: Overwriting existing fragment BED %s. Continue? [y/N] " %
            output_fp)
        if not overwrite.lower() == 'y':
            print("Did not overwrite existing fragment BED.")
            return
        os.remove(output_fp)

    print("Digesting")
    if "fasta" in genome or "fa" in genome:
        genome = Bio.SeqIO.parse(open(genome, "rU"), format='fasta')
    else:
        genome = Bio.SeqIO.parse(open(genome, "rU"), format='genbank')

    with open(output_fp, 'w') as bedfile:
        output = csv.writer(bedfile, delimiter="\t")
        for chromosome in genome:
            print(chromosome.id, chromosome.name)
            # Digest the sequence data and return the cut points
            fragment_num = 0
            enzyme = RestrictionBatch([restriction_enzyme])
            for enzyme, cutpoints in enzyme.search(
                    chromosome.seq, linear=linear).items():
                cutpoints.insert(0, 0)
                # Covert cut points to fragments
                for index, point in enumerate(cutpoints):
                    # Adjust for start and end of sequence and the offset for the cutpoint
                    # ATTENTION: This will only work for MspI I will have to spend some time to make it compatible with
                    #		   any restriction enzyme such as ones with blunt ends
                    #startpoint = 0 if cutpoints[index] - 1 < 0 else cutpoints[index] - 1
                    startpoint = 0 if cutpoints[index] - \
                        1 < 0 else cutpoints[index]
                    endpoint = len(chromosome.seq) if index + \
                        1 >= len(cutpoints) else cutpoints[index + 1] - 1

                    accession = chromosome.id
                    version = ''
                    if "." in chromosome.id:
                        accession, version = chromosome.id.split(".")
                    if not accession.startswith("chr"):
                        accession = "chr" + accession
                    output.writerow(
                        [accession, startpoint, endpoint, fragment_num])
                    fragment_num += 1
    if not do_not_index:
        build_fragment_index(output_fp, output_db)


def build_fragment_index(fragment_fp, output_db):
    if not output_db:
        if not fragment_fp.rfind('.') == -1:
            output_db = fragment_fp[:fragment_fp.rfind('.')] + ".db"
        else:
            output_db = fragment_fp + ".db"
    if os.path.isfile(output_db):
        overwrite = input(
            "WARNING: Overwriting existing fragment database %s. Continue? [y/N] " %
            output_db)
        if not overwrite.lower() == 'y':
            print("Did not overwrite existing fragment database.")
            return
        os.remove(output_db)

    if not os.path.isdir(os.path.dirname(output_db)):
        os.path.makedirs(os.path.dirname(output_db))

    fragment_index_db = sqlite3.connect(output_db)
    fragment_index = fragment_index_db.cursor()
    fragment_index.execute(
        "CREATE TABLE fragments (chr text, start integer, end integer, fragment integer)")
    fragment_index.execute("CREATE INDEX f_index ON fragments (chr,fragment)")

    with open(fragment_fp, 'r') as fragments_bed:
        for line in fragments_bed:
            fragment = line.strip().split('\t')
            fragment_index.execute("INSERT INTO fragments VALUES (?,?,?,?)", [
                                   fragment[0][fragment[0].find("chr") + 3:], int(fragment[1]), int(fragment[2]), fragment[3]])
    fragment_index_db.commit()


def build_gene_index(
        gene_files,
        output_bed,
        output_db,
        config,
        symbol_col=27,
        chr_col=29,
        start_col=30,
        end_col=31,
        p_thresh_col=None,
        no_header=False,
        do_not_tidy_up=False):
    genes = {}

    # Edited to cater for when only the .gtf and .conf arguments are provided.
    if symbol_col:
        symbol_col -= 1
    if chr_col:
        chr_col -= 1
    if start_col:
        start_col -= 1
    if end_col:
        end_col -= 1
    if p_thresh_col:
        p_thresh_col -= 1

    if not output_bed:
        output_bed = os.path.join(
            config.get(
                "Defaults",
                "LIB_DIR"),
            "gene_reference.bed")

    if not output_db:  # corrected output_bed to output_db
        output_db = os.path.join(
            config.get(
                "Defaults",
                "LIB_DIR"),
            "gene_reference.db")

    append_bed = False
    overwrite_bed = True
    if os.path.isfile(output_bed):
        upsert = input(
            "WARNING: Appending input to existing BED file (%s). Continue? [y/N] " %
            output_bed)
        if not upsert.lower() == 'y':
            print("Did not append to existing gene database.")
            overwrite_bed = False
        else:
            append_bed = True

    if (overwrite_bed or append_bed) and not os.path.isdir(
            os.path.dirname(output_bed)):
        os.makedirs(os.path.dirname(output_bed))

    upsert_db = True

    if os.path.isfile(output_db):
        upsert = input(
            "WARNING: Upserting input to existing gene database (%s). Continue? [y/N] " %
            output_db)
        if not upsert.lower() == 'y':
            print("Did not write to existing gene database.")
            upsert_db = False

    if upsert_db and not os.path.isdir(os.path.dirname(output_db)):
        os.makedirs(os.path.dirname(output_db))

    if (not (overwrite_bed or append_bed)) and not upsert_db:
        print("No action performed; exiting.")
        return

    if upsert_db:
        gene_index_db = sqlite3.connect(output_db)
        gene_index = gene_index_db.cursor()
        if p_thresh_col:
            gene_index.execute(
                "CREATE TABLE IF NOT EXISTS genes (symbol text, chr text, start integer, end integer, double p_thresh)")
        else:
            gene_index.execute(
                "CREATE TABLE IF NOT EXISTS genes (symbol text, chr text, start integer, end integer)")
        gene_index.execute(
            "CREATE INDEX IF NOT EXISTS g_index ON genes (symbol)")

    for gene_file in gene_files:
        # Do line count for progress meter
        lines = 0
        with open(gene_file, 'r') as genefile:
            print("Determining table size...")
            for i in genefile:
                lines += 1
            lines = lines // 100 * 100  # Get an approximation
            do_linecount = not lines == 0

        with open(gene_file, 'r') as genefile:
            # Determine if input file is a GTEx-supplied gene reference
            is_gtex_file = gene_file.endswith(".gtf")
            # If so, process headers accordingly
            if is_gtex_file:
                line = genefile.readline()
                while line.startswith("##"):
                    line = genefile.readline()
            # Otherwise, process headers according to script options
            else:
                if not no_header:
                    genefile.readline()
                line = genefile.readline()
            i = 0
            # For each line
            while line:
                if do_linecount:
                    if i % (lines / 100) == 0:
                        print("\tProcessed %d%%..." % ((float(i) / float(lines)) * 100))
                # If the line is a GTEx file, extract information accordingly
                if is_gtex_file:
                    entry = line.strip().split('\t')
                    if entry[2] == "gene":
                        gene_stats = entry[8].split(';')
                        gene_symbol = gene_stats[4].strip().split(' ')[
                            1].strip('"')
                        gene_chr = entry[0]
                        gene_start = int(entry[3])
                        gene_end = int(entry[4])
                    else:
                        i += 1
                        line = genefile.readline()
                        continue  # Skip if the entry is not for a canonical gene
                # Otherwise, extract by column
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
                # Enter into index, regardless of input file type
                if not gene_symbol in genes:
                    genes[gene_symbol] = {
                        "chr": gene_chr, "start": gene_start, "end": gene_end}
                    if p_thresh_col:
                        genes[gene_symbol]["p_thresh"] = gene_p_thresh
                else:
                    curr_length = genes[gene_symbol]["end"] - \
                        genes[gene_symbol]["start"]
                    if curr_length < abs(gene_end - gene_start):
                        genes[gene_symbol]["start"] = gene_start
                        genes[gene_symbol]["end"] = gene_end
                line = genefile.readline()
                i += 1

    bed_out = None
    if overwrite_bed:
        bed_out = open(output_bed, 'w')
    elif append_bed:
        bed_out = open(output_bed, 'a')
    for gene in genes:
        if bed_out:
            bed_out.write(
                'chr%s\t%s\t%s\t%s\n' %
                (genes[gene]["chr"],
                 genes[gene]["start"],
                    genes[gene]["end"],
                    gene))
        if upsert_db:
            if p_thresh_col:
                gene_index.execute(
                    "INSERT INTO genes VALUES (?,?,?,?,?)", [
                        gene, genes[gene]["chr"], genes[gene]["start"], genes[gene]["end"], genes[gene]["p_thresh"]])
            else:
                gene_index.execute(
                    "INSERT INTO genes VALUES (?,?,?,?)", [
                        gene, genes[gene]["chr"], genes[gene]["start"], genes[gene]["end"]])
    if bed_out:
        bed_out.close()
    if upsert_db:
        gene_index_db.commit()
        gene_index.close()

    if not do_not_tidy_up:
        print("Tidying up...")
        for gene_file in gene_files:
            os.remove(gene_file)


def build_eqtl_index(
        table_fp,
        output_fp=None,
        snp_col=23,
        gene_symbol_col=27,
        gene_chr_col=29,
        gene_start_col=30,
        gene_stop_col=31,
        p_val_col=6,
        effect_size_col=3,
        do_not_tidy_up=False):
    if not output_fp:
        if not table_fp.rfind('.') == -1:
            output_fp = table_fp[:table_fp.rfind('.')] + ".db"
        else:
            output_fp = table_fp + ".db"

    if not os.path.isdir(os.path.dirname(output_fp)):
        os.makedirs(os.path.dirname(output_fp))

    if os.path.isfile(output_fp):
        upsert = input(
            "WARNING: Upserting input to existing eQTL database %s. Continue? [y/N] " %
            output_fp)
        if not upsert.lower() == 'y':
            print("Did not write to existing eQTL database.")
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
    table_index.execute(
        "CREATE TABLE IF NOT EXISTS eqtls (rsID text, gene_name text, gene_chr text, gene_start integer, gene_end integer, pvalue real, effect_size real)")
    table_index.execute("CREATE INDEX IF NOT EXISTS id ON eqtls (rsID)")
    # Do line count for progress meter
    do_linecount = True
    print("Determining table size...")
    with open(table_fp, 'r') as eqtl_table:
        lines = 0
        for i in eqtl_table:
            lines += 1
    lines = lines // 100 * 100  # Get an approximation
    do_linecount = not lines == 0

    with open(table_fp, 'r') as eqtl_table:
        for i, line in enumerate(eqtl_table):
            if do_linecount:
                if i % (lines / 100) == 0:
                    print("\tProcessed %d%%..." % ((float(i) / float(lines)) * 100))
            if i == 0:
                continue
            eqtl = line.strip().split('\t')
            table_index.execute(
                "INSERT INTO eqtls VALUES (?,?,?,?,?,?,?)",
                [
                    eqtl[snp_col],
                    eqtl[gene_symbol_col],
                    eqtl[gene_chr_col],
                    eqtl[gene_start_col],
                    eqtl[gene_stop_col],
                    eqtl[p_val_col],
                    eqtl[effect_size_col]])
    table_index_db.commit()
    print("Done indexing eQTL table.")
    if not do_not_tidy_up:
        print("Tidying up...")
        os.remove(table_fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--inputs", nargs='+', required=True,
        help="The the dbSNP IDs or loci of SNPs of interest in the format " +\
                        "\"chr<x>:<locus>\"")
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used for this " +\
                        "hiCquery run (default: conf.py)")
    parser.add_argument("-n", "--include_cell_lines", nargs='+',
                        help="Space-separated list of cell lines to include " +\
                        "(others will be ignored). NOTE: Mutually exclusive " +\
                        "with EXCLUDE_CELL_LINES.")
    parser.add_argument("-x", "--exclude_cell_lines", nargs='+',
                        help="Space-separated list of cell lines to exclude " +\
                        "(others will be included). NOTE: Mutually exclusive " +\
                        "with INCLUDE_CELL_LINES.")
    parser.add_argument("-o", "--output_dir", default="codes3d_output",
                        help="The directory in which to output results " +\
                        "(\"hiCquery_output\" by default).")
    parser.add_argument("-l", "--local_databases_only", action="store_true",
                        default=False, help="Consider only local databases. " +\
                        "Will only include cis-eQTLs if using downloadable " +\
                        "GTEx dataset.")
    parser.add_argument("-s", "--suppress_intermediate_files", action="store_true",
                        default=False, help="Do not produce intermediate " +\
                        "files. These can be used to run the pipeline from " +\
                        "an intermediate stage in the event of interruption.")
    parser.add_argument("-p", "--num_processes", type=int, default=1,
                        help="Desired number of processes for multithreading " +\
                        "(default: 1).")
    parser.add_argument("-f", "--fdr_threshold", type=float, default=0.05,
                        help="The FDR threshold to consider an eQTL " +\
                        "statistically significant (default: 0.05).")
    parser.add_argument("-r", "--restriction_enzymes", nargs='+',
                        help="Space-separated list of  " +\
                        "restriction enzymes used in HiC data.")
    parser.add_argument("-b","--buffer_size_in",type=int,default=1048576,
                        help="Buffer size applied to file input during " +\
                        "compilation (default: 1048576).")
    parser.add_argument("-d","--buffer_size_out",type=int,default=1048576,
                        help="Buffer size applied to file output during " +\
                        "compilation (default: 1048576).")
    parser.add_argument("-t", "--num_processes_summary", type=int, 
                        default=min(psutil.cpu_count(), 32),
                        help="The number of processes for compilation of " +\
                        "the results (default: %s)." %
                        str(min(psutil.cpu_count(), 32)))
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    lib_dir = os.path.join(os.path.dirname(__file__),config.get("Defaults", "LIB_DIR"))
    snp_database_fp = os.path.join(os.path.dirname(__file__),
                                   config.get("Defaults", "SNP_DATABASE_FP"))
    hic_data_dir = os.path.join(os.path.dirname(__file__),
                                config.get("Defaults", "HIC_DATA_DIR"))
    fragment_bed_fp = os.path.join(os.path.dirname(__file__),
                                   config.get("Defaults", "FRAGMENT_BED_FP"))
    fragment_database_fp = os.path.join(os.path.dirname(__file__),
                                        config.get("Defaults", "FRAGMENT_DATABASE_FP"))
    gene_bed_fp = os.path.join(os.path.dirname(__file__),
                               config.get("Defaults", "GENE_BED_FP"))
    gene_database_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults", "GENE_DATABASE_FP"))
    eqtl_data_dir = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "EQTL_DATA_DIR"))
    expression_table_fp = os.path.join(os.path.dirname(__file__),
                                       config.get("Defaults", "EXPRESSION_TABLE_FP"))
    gene_dict_fp = os.path.join(os.path.dirname(__file__),
                                       config.get("Defaults", "GENE_DICT_FP"))
    snp_dict_fp = os.path.join(os.path.dirname(__file__),
                                       config.get("Defaults", "SNP_DICT_FP"))
    GTEX_CERT = os.path.join(os.path.dirname(__file__),
                             config.get("Defaults", "GTEX_CERT"))
    HIC_RESTRICTION_ENZYMES = [e.strip() for e in \
                                config.get("Defaults", "HIC_RESTRICTION_ENZYMES").split(',')]
    restriction_enzymes, include_cell_lines, exclude_cell_lines =\
        parse_parameters(args.restriction_enzymes, \
                         args.include_cell_lines, \
                         args.exclude_cell_lines)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    snps = process_inputs(
        args.inputs, snp_database_fp, lib_dir,
        restriction_enzymes, args.output_dir,
        args.suppress_intermediate_files)
    interactions = find_interactions(
        snps, lib_dir, include_cell_lines,
        exclude_cell_lines, args.output_dir,
        args.suppress_intermediate_files)
    genes = find_genes(interactions, lib_dir, gene_bed_fp,
                       args.output_dir, args.suppress_intermediate_files)
    num_sig, p_values = find_eqtls(
        snps, genes, eqtl_data_dir, gene_database_fp,
        args.fdr_threshold, args.local_databases_only,
        args.num_processes, args.output_dir, gene_dict_fp, snp_dict_fp,
        suppress_intermediate_files=args.suppress_intermediate_files)
    produce_summary(
         p_values, snps, genes, gene_database_fp, expression_table_fp,
         args.fdr_threshold, args.output_dir, args.buffer_size_in,
         args.buffer_size_out, args.num_processes_summary)
    #produce_overview(genes,eqtls,num_sig,args.output_dir)
    #pathways = retrieve_pathways(eqtls,args.fdr_threshold,args.num_processes,args.output_dir)

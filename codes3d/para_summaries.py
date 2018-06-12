#!usr/bin/python

import csv
import random
import os
import sqlite3
import multiprocessing
import shutil
import codes3d
import configparser

"""
Date: 23052017
Purpose: To multiprocess the process of summaryy files for negative controls
"""

def initialise_traits(traits_file, snps_dir):
    trait_snps = {}
    traits_file = open(traits_file, 'r')
    traits = traits_file.readlines()
    traits = [trait.strip() for trait in traits]

    if os.path.isdir(snps_dir):
        for trait in traits:
            tfile = open(os.path.join(snps_dir, trait) + '.txt', 'r')
            snps = tfile.readlines()
            snps = [snp.strip() for snp in snps]
            trait_snps[trait] = snps

    return trait_snps

def generate_replicates(trait_snps, random_daSNPs_dir):
    snps = []
    snp_array = {}
    for trait in trait_snps:
        snps += trait_snps[trait]

    for i in range(0, 1000):
        snp_list = []
        for r in range(0, 483):
            rand_snp = random.choice(snps)
            while rand_snp in snp_list:
                rand_snp = random.choice(snps)
            snp_list.append(rand_snp)
        snp_array[i] = snp_list

    if not os.path.isdir(random_daSNPs_dir):
        os.mkdir(random_daSNPs_dir)
    for rep in snp_array:
        sfile = open(os.path.join(random_daSNPs_dir, str(rep)) + '.txt', 'w')
        for snp in snp_array[rep]:
            sfile.write(snp + '\n')
        sfile.close()

def retrieve_eqtls(rep):
    """
    random_daSNPs_dir, eqtlDB = params
    conn = sqlite3.connect(eqtlDB)
    conn.text_factory = str
    cur = conn.cursor()
    
    _,_,replicates = next(os.walk(random_daSNPs_dir), (None, None, []))

    all_eqtls = {}
    
    print('Retrieving eqtls for...')
    for rep in replicates:
        rep_path = os.path.join(random_daSNPs_dir, rep)
        rfile = open(rep_path, 'r')
        snps = rfile.readlines()
        rfile.close()
        all_eqtls[rep] = {}
        snps = [snp.strip() for snp in snps]
        non_eqtls = []
        rep_dir = (rep_path[:len(rep_path)-4])
        print('\t replicate set ' + rep[:len(rep)-4])        
        if not os.path.isdir(rep_dir):
            os.mkdir(rep_dir)
        for snp in snps:
            cur.execute("SELECT * FROM eqtls WHERE snp = ?;", (snp,))
            snp_eqtls = cur.fetchall()
            if snp_eqtls:
                all_eqtls[rep][snp] = snp_eqtls
                print('\t\t' + str(snps.index(snp)+1) + '\t' + snp + \
                          ' (' + str(len(snp_eqtls)) + ' eqtls)')
            else:
                print('\t\t' +  str(snps.index(snp)+1) + '\t' + snp + \
                          ' not found in eqtl database')
                non_eqtls.append(snp)
        eqtl_file = open(rep_dir + '/eqtls.txt', 'w')
        ewriter = csv.writer(eqtl_file, delimiter = '\t')
        print('\t\t Writing eqtl files')
        for snp in all_eqtls[rep]:
            ewriter.writerows(all_eqtls[rep][snp])
        eqtl_file.close()
        non_eqtl_file = open(rep_dir + '/non-eqtls.txt', 'w')
        for snp in non_eqtls:
            non_eqtl_file.write(snp + '\n')
        non_eqtl_file.close()
        shutil.move(rep_path, rep_dir + rep)
    """
    #random_daSNPs_dir, rep = params
    rep_path = os.path.join(random_daSNPs_dir, rep)

    conn = sqlite3.connect(eqtlDB)
    conn.text_factory = str
    cur = conn.cursor()
    
    all_eqtls = {}
    if os.path.isfile(rep_path):
        rfile = open(rep_path, 'r')
        snps = rfile.readlines()
        rfile.close()
        all_eqtls[rep] = {}
        snps = [snp.strip() for snp in snps]
        non_eqtls = []
        rep_dir = (rep_path[:len(rep_path)-4])
        print('\t replicate set ' + rep[:len(rep)-4])
        if not os.path.isdir(rep_dir):
            os.mkdir(rep_dir)
        for snp in snps:
            cur.execute("SELECT * FROM eqtls WHERE snp = ?;", (snp,))
            snp_eqtls = cur.fetchall()
            if snp_eqtls:
                all_eqtls[rep][snp] = snp_eqtls
                print('\t\t' + 'rep' + rep[:len(rep)-4] + ':' + \
                          str(snps.index(snp)+1) + '\t' + snp + \
                          ' (' + str(len(snp_eqtls)) + ' eqtls)')
            else:
                print('\t\t' +  'rep' + rep[:len(rep)-4] + ':' + \
                          str(snps.index(snp)+1) + '\t' + snp + \
                          ' not found in eqtl database')
                non_eqtls.append(snp)
        eqtl_file = open(rep_dir + '/eqtls.txt', 'w')
        ewriter = csv.writer(eqtl_file, delimiter = '\t')
        print('\t\t Writing eqtl files')
        for snp in all_eqtls[rep]:
            ewriter.writerows(all_eqtls[rep][snp])
        eqtl_file.close()
        non_eqtl_file = open(rep_dir + '/non-eqtls.txt', 'w')
        for snp in non_eqtls:
            non_eqtl_file.write(snp + '\n')
        non_eqtl_file.close()
    
        shutil.move(rep_path, rep_dir + '/' + rep)
    else:
        print(rep, ': File not found :)')


def correct(replicates):
    print(len(replicates))
    i = 20
    for rep in replicates:
        #shutil.move(os.path.join(random_daSNPs_dir, rep), \
        #                os.path.join(random_daSNPs_dir, (str(i)+'.txt')))
        i += 1
        print(rep)

def file_len(fname):
    i = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def col_len(fname):
    f = open(fname, 'r')
    creader = csv.reader(f, delimiter = '\t')
    for row in creader:
        if len(row) < 11:
            return(false)

def prod_summary(repdir):
    daSNPs_dir = '/mnt/3dgenome/projects/tfad334/hvn/analysis/'
    daSNPs_dir += 'control/random_daSNPs'
    eqtl_file = os.path.join(daSNPs_dir, repdir) + '/eqtls.txt'
    if file_len(eqtl_file): #Handle traits with no eQTLs
        config_file = 'docs/codes3d.conf'
        output_dir = os.path.join(daSNPs_dir, repdir)
        config = configparser.ConfigParser()
        config.read(config_file)
        expression_table_fp = config.get("Defaults","EXPRESSION_TABLE_FP")
        print(repdir)
        eqtls,num_sig = codes3d.parse_eqtls_files(eqtl_file, 0.05)
        codes3d.produce_summary(eqtls,expression_table_fp,output_dir)
    else:
        return(eqtl_file)
if __name__ == '__main__':
    daSNPs_dir = '/mnt/3dgenome/projects/tfad334/hvn/analysis/'
    daSNPs_dir += 'control/random_daSNPs'
    pool = multiprocessing.Pool(processes=20)
    _,repdirs,_ = next(os.walk(daSNPs_dir), (None, [],None))
    repdirs = [repdir for repdir in repdirs 
               if (os.path.isfile(os.path.join(daSNPs_dir,repdir) + '/eqtls.txt')
               and os.path.isfile(os.path.join(daSNPs_dir, repdir) + 
                                      '/summary.txt'))
               ]
    pool.map(prod_summary, repdirs)
    #problem_traits = []
    #for repdir in repdirs:
    #    try:
    #        prod_summary(repdir)
    #    except IndexError:
    #        problem_traits.append(repdir)
    #print(problem_traits)

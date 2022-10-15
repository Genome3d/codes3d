#!/usr/bin/env python

import pandas as pd
import os
import sys
import sqlite3
import argparse
from sqlalchemy import create_engine
import csv


def read_promoter_file(promoters):
    col = ['frag_chr', 'frag_start', 'frag_end', 'frag_id', 'gene', 'project']
    gene_prom = pd.read_csv(promoters, sep = "\t", names = col, 
            usecols = ['frag_id', 'gene', 'project'])
    return gene_prom

def read_genebed(genebed):
    cols = ['chr', 'start', 'end', 'gene', 'gencode_id']
    gene_bed = pd.read_csv(genebed, sep = "\t", names = cols)
    return gene_bed

def build_gene_table(gene_prom, gene_bed, table, db_url, outfile):
    merge_tables = pd.merge(gene_prom, gene_bed, how='inner', on='gene', sort=False)
    gene_file= merge_tables[['chr', 'start', 'end', 'gene', 'gencode_id', 
        'frag_id', 'project']].drop_duplicates()
    db = create_engine(db_url, echo=False)
    gene_file.to_csv(outfile, sep='\t', header=True, index=False)
    gene_file.to_sql(table, con=db, if_exists='replace', index=False)
    print("Creating index...")
    db.execute("CREATE INDEX idx_{}_frag_id ON {}(frag_id)".format(table, table))
    db.execute('''CREATE INDEX idx_{}_fragid_gencodeid on 
            {}(frag_id,gencode_id)'''.format(table, table))
    db.execute('''CREATE INDEX idx_{}_fragid_gencodeid_project on
            {}(frag_id,gencode_id,project)'''.format(table, table))
    db.execute('''CREATE INDEX idx_{}_gencodeid_project on
            {}(gencode_id,project)'''.format(table, table)) 
    db.execute('''CREATE INDEX idx_{}_gencode_id on 
            {}(gencode_id)'''.format(table, table))
    print("Done.")
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 
            'Build gene tables for promoter capture Hi-C data analysis')
    parser.add_argument('-p', dest = 'pro_file', nargs='+', required = True,
            help = '''Provide a path to the promoter fragment file (no header).
            Expected columns in the file: frag_chr, frag_start, frag_end, frag_id, gene_name,
            project. Multiple files should be speparated by a space''') 
    parser.add_argument('-g', dest = 'gene_bed', required = True,
            help='''Path to the gene bed file with extension .bed.
            e.g. /mnt/projects/codes3d/lib/reference_files/genes/gene_reference.bed''')
    parser.add_argument("-t", dest = 'table', required = True, 
            help='Name of table e.g. gene_lookup_pchic_hindiii.')
    parser.add_argument("-u", dest = 'db_url', required = True,
            help="URL of database e.g. postgresql://user:password@hostname/database")
    parser.add_argument("-o", dest = 'outfile', required = True,
            help="absolute path to the output file")
    args = parser.parse_args()

    promoters = []
    if len(args.pro_file) > 1:
        for pfile in args.pro_file:
            p_df = read_promoter_file(pfile)
            promoters.append(p_df)
        promoters = pd.concat(promoters)
    else:
        promoters = read_promoter_file(args.pro_file[0])
    
    genebed = read_genebed(args.gene_bed)

    build_gene_table(promoters, 
            genebed, 
            args.table, 
            args.db_url, 
            args.outfile)
   


#! /usr/bin/env python

import torch
from tensorqtl import *
from tensorqtl import genotypeio, trans
from tensorqtl import tensorqtl
import argparse
import pandas as pd
from sqlalchemy import create_engine
import time
from tqdm import tqdm
sys.path.insert(
    1, os.path.join(os.path.abspath(os.path.dirname(__file__)), 'tensorqtl'))

'''Create genotype table for eQTL mapping by tensorqtl
'''


def create_genotype_table(genotype_df, db_url):
    table = 'genotype'
    desc = 'Building genotype table...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='variants', desc=desc, disable=False,
             bar_format=bar_format)
    chunksize = 50000
    chunks = [genotype_df[i:i+chunksize]
              for i in range(0, len(genotype_df), chunksize)]
    db = create_engine(db_url, echo=False)
    for df in chunks:
        if_exists = 'replace' if t.total == 0 else 'append'
        df.to_sql(table, con=db, if_exists=if_exists)
        t.total += len(df)
        t.update(len(df))
    t.close()
    # db.execute('''CREATE INDEX idx_{}_snp ON {}
    # (snp)'''.format(table, table))


def read_genotype(plink_prefix):
    ''' Read WGS VCF '''
    pr = genotypeio.PlinkReader(plink_prefix)
    genotype_df = pd.DataFrame(pr.load_genotypes(),
                               #index=pr.bim['snp'],
                               columns=pr.fam['iid'])
    return genotype_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-p', '--plink-prefix', required=True,
        help='''Filepath to plink generated files of genotypes but without the 
        extensions''')
    parser.add_argument(
        '-u', '--db-url', required=True,
        help='URL of database e.g postgresql://user:password@hostname/database')
    args = parser.parse_args()
    start_time = time.time()
    genotype_df = read_genotype(args.plink_prefix)
    print(genotype_df)
    create_genotype_table(genotype_df, args.db_url)
    print('Time elapsed: {:.2f} mins'.format(
        (time.time() - start_time)/60))

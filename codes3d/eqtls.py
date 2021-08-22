#! /usr/bin/env python

from tensorqtl import trans
import tensorqtl
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import os
import sys
import pandas as pd
import numpy as np
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
import time
import multiprocessing
import tqdm
from itertools import repeat
#import sqlite3
import statsmodels.stats.multitest as multitest
#from sqlalchemy_utils import database_exists
import argparse
import configparser
import codes3d

sys.path.insert(
    1, os.path.join(os.path.abspath(os.path.dirname(__file__)), 'tensorqtl'))


def create_db(tissue, logger):
    conn = None
    database = 'eqtls_gtex_{}'.format(tissue.lower())
    try:
        conn = psycopg2.connect(
            host='127.0.0.1',
            # database=database,
            user='codes3d',
            password='codes3d'
        )
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        cur = conn.cursor()
        cur.execute(
            '''DROP DATABASE IF EXISTS {}; '''.format(database)
        )
        cur.execute(
            '''CREATE DATABASE {} TABLESPACE {}; '''.format(
                database, 'tblspace_codes3d_commons')
        )
        cur.close()
        logger.write('{} DATABASE created.'.format(database))
    except (Exception, psycopg2.DatabaseError) as error:
        logger.write(error)
    finally:
        if conn is not None:
            conn.close()
            create_eqtls_table(tissue)


def create_eqtls_table(tissue, logger):
    conn = None
    database = 'eqtls_gtex_{}'.format(tissue.lower())
    try:
        conn = psycopg2.connect(
            host='127.0.0.1',
            database=database,
            user='codes3d',
            password='codes3d'
        )
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        cur = conn.cursor()
        cur.execute('''DROP TABLE IF EXISTS {}'''.format(database))
        cur.execute(
            '''CREATE TABLE eqtls(
            variant_id TEXT,
            phenotype_id TEXT,
            pval DOUBLE PRECISION,
            maf DOUBLE PRECISION
            );
            '''
        )
        cur.execute(
            '''CREATE INDEX idx_eqtls ON eqtls(variant_id, phenotype_id); '''
        )
        cur.close()
        logger.write('eqtls TABLE created in {}.'.format(database))
    except (Exception, psycopg2.DatabaseError) as error:
        logger.write(error)
    finally:
        if conn is not None:
            conn.close()


def migrate_to_postgres(fp, postgres_url):
    for tissue_fp in os.listdir(fp):
        tissue = tissue_fp.split('.')[0]
        if not tissue.startswith('A'):
            continue
        start_time = time.time()
        create_db(tissue)
        postgres_db = create_engine(
            postgres_url.format('eqtls_gtex_' + tissue.lower()),
            echo=False)
        sqlite_db = create_engine(
            'sqlite:////{}'.format(os.path.join(fp, tissue_fp)), echo=False)
        sql = '''SELECT * FROM eqtls'''
        for df in pd.read_sql_query(sql, con=sqlite_db, chunksize=100000):
            df.to_sql('eqtls', con=postgres_db,
                      if_exists='append', index=False)
        logger.write('Time elapsed: {} mins.\n'.foramt(
            (time.time()-start_time)/60))


def map_tissue_eqtls(
        tissue,
        pairs_df,
        genotype_df,
        variant_df,
        eqtls,
        covariates_dir,
        expression_dir,
        pval_threshold,
        maf_threshold,
        eqtl_project
):
    covariates_fp = ''
    phenotype_fp = ''
    if eqtl_project.lower() == 'gtex': # TODO rename files to be consistent.
        covariates_fp = os.path.join(
            covariates_dir, tissue + '.v8.covariates.txt')
        phenotype_fp = os.path.join(
            expression_dir, tissue + '.v8.normalized_expression.bed.gz')
    else:
        covariates_fp = os.path.join(
            covariates_dir, tissue + '.covariates.txt')
        phenotype_fp = os.path.join(
            expression_dir, tissue + '.normalized_expression.bed.gz')
    if not (os.path.exists(covariates_fp) and os.path.exists(phenotype_fp)):
        return
    covariates_df = pd.read_csv(covariates_fp, sep='\t', index_col=0).T
    phenotype_df, pos_df = tensorqtl.read_phenotype_bed(phenotype_fp)
    phenotype_df = phenotype_df[
        phenotype_df.index.isin(pairs_df['pid'])]
    eqtl_df = trans.map_trans(
        genotype_df,
        phenotype_df,
        covariates_df,
        return_sparse=True,
        pval_threshold=float(pval_threshold),
        maf_threshold=float(maf_threshold),
        batch_size=20000,
        verbose=False)
    eqtl_df['tissue'] = tissue
    eqtls.append(eqtl_df[~((eqtl_df['variant_id'].isnull()) |
                           (eqtl_df['phenotype_id'].isnull()))])

def fetch_genotypes_deprecated(db, df):
    genotypes = []
    sql = '''SELECT * FROM genotype WHERE snp = '{}' '''
    with db.connect() as con:
        for snp in df.index.values:
            genotypes.append(pd.read_sql(
                sql.format(snp), con=con))
    if len(genotypes) > 0:
        return pd.concat(genotypes)
    return pd.DataFrame()


def fetch_genotypes(db, snp):
    df = pd.DataFrame()
    with db.connect() as con:
        df = pd.read_sql(
            '''SELECT * FROM genotype WHERE snp = '{}' '''.format(snp)
            , con=con)
    if not df.empty:
        return df.iloc[0]
    else:
        return pd.DataFrame({'snp': [snp]}, columns=df.columns).iloc[0]

def chunk_list(it, size):
    for i in range(0, len(it), size):
        yield it[i:i+size]
        
def calc_chunksize(loci, distance):
    if len(loci) == 1:
        return 1
    diff = abs(max(loci) - min(loci))
    while diff > distance:
        chunks = [chunk for chunk in chunk_list(loci, int((len(loci)/2)+1))]
        if(len(chunks)==1):
            return(diff)
        if abs(max(chunks[0])-min(chunks[0])) > abs(max(chunks[1])-min(chunks[1])):
            loci = sorted(chunks[0])
        else:
            loci = sorted(chunks[1])
        diff = abs(max(loci)-min(loci))
    return diff
            
def fetch_genotypes_beta(db, df):
    genotypes = []
    chroms = df.groupby(['chrom']).agg({'pos': [np.min, np.max]})
    chroms.columns = chroms.columns.droplevel()
    chroms['distance'] = chroms['amax'] - chroms['amin']
    chunksize = calc_chunksize(
        df[df['chrom']==chroms['distance'].idxmax()].sort_index().reset_index(
            drop=True)['pos'].astype(int).tolist(),
        20000 # Bigger numbers can cause out shared memory errors
    )
    chunks = [df[i:i+chunksize] for i in range(0, len(df), chunksize)]
    desc = '  * Loading genotypes'
    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'

    pbar = tqdm.tqdm(total=len(df), desc=desc, bar_format=bar_format,
                     unit='variants', ncols=80)
    sql = '''SELECT * FROM genotype WHERE snp >= '{}' AND snp <= '{}' '''
    with db as con:
        for chunk in chunks:
            res = pd.read_sql(
                sql.format(chunk.index.min(), chunk.index.max()), con=con)
            genotypes.append(res[res['snp'].isin(chunk.index)])
            pbar.update(len(chunk))
    pbar.close()
    if len(genotypes) > 0:
        return pd.concat(genotypes)
    return pd.DataFrame()


def fetch_genotypes_raw(db, df):
    genotypes = []
    desc = '  * Loading genotypes'
    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
    with db.connect() as con:
        for idx, row in tqdm.tqdm(
                df.iterrows(), total=len(df), desc=desc, bar_format=bar_format,
                unit='variants', ncols=80):
            res = pd.read_sql(
                '''SELECT * FROM genotype WHERE snp = '{}' '''.format(
                    row.name), con=con)
            if not res.empty:
                genotypes.append(res)
    
    if len(genotypes) > 0:
        return pd.concat(genotypes)
    else:
        return pd.DataFrame()

def map_eqtls(
        gene_df,
        tissues,
        num_processes,
        db,
        covariates_dir,
        expression_dir,
        pval_threshold,
        maf_threshold,
        fdr_threshold,
        logger
):
    start_time = time.time()
    logger.write("Identifying eQTLs...")
    eqtl_project = tissues['project'].tolist()[0]
    tissues = tissues['name'].tolist()
    pairs_df = gene_df[
        ['variant_id', 'gencode_id', 'chrom', 'locus']
    ].drop_duplicates()
    pairs_df.columns = ['sid', 'pid', 'sid_chr', 'sid_pos']
    variant_df = pairs_df[['sid', 'sid_chr', 'sid_pos']].drop_duplicates()
    variant_df = variant_df.rename(columns={
        'sid': 'snp',
        'sid_chr': 'chrom',
        'sid_pos': 'pos'})
    variant_df = variant_df.set_index('snp')
    genotype_df = fetch_genotypes_raw(db, variant_df)
    if genotype_df.empty:
        logger.write('Variants were not found in database.')
        return
    genotype_df.set_index('snp', inplace=True)
    genotype_df.columns.name = 'iid'
    manager = multiprocessing.Manager()
    eqtls = manager.list()
    desc = '  * Mapping eQTLs'
    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
    if len(tissues) > 1:
        with multiprocessing.Pool(num_processes) as pool:
            for _ in tqdm.tqdm(
                    pool.istarmap(
                        map_tissue_eqtls,
                        zip(tissues,
                            repeat(pairs_df),
                            repeat(genotype_df),
                            repeat(variant_df),
                            repeat(eqtls),
                            repeat(covariates_dir),
                            repeat(expression_dir),
                            repeat(pval_threshold),
                            repeat(maf_threshold),
                            repeat(eqtl_project))
                    ),
                    total=len(tissues), desc=desc, bar_format=bar_format,
                    unit='tissues', ncols=80):
                pass
    else:
        batch_size = int(len(pairs_df)/num_processes)
        if batch_size < 1:
            batch_size = 1
        pairs_batches = [pairs_df[i:i+batch_size]
                         for i in range(0, len(pairs_df), batch_size)]        
        with multiprocessing.Pool(num_processes) as pool:
            for _ in tqdm.tqdm(
                    pool.istarmap(
                        map_tissue_eqtls,
                        zip(repeat(tissues[0]),
                            pairs_batches,
                            repeat(genotype_df),
                            repeat(variant_df),
                            repeat(eqtls),
                            repeat(covariates_dir),
                            repeat(expression_dir),
                            repeat(pval_threshold),
                            repeat(maf_threshold),
                            repeat(eqtl_project))
                    ),
                    total=len(pairs_batches), desc=desc, bar_format=bar_format,
                    unit='batches', ncols=80):
                pass
    eqtl_df = pd.DataFrame()
    if len(eqtls) > 0:
        eqtl_df = pd.concat(eqtls)
    if eqtl_df.empty:
        return eqtl_df
    eqtl_df = eqtl_df.merge(
        pairs_df, how='inner',
        right_on=['sid', 'pid'],
        left_on=['variant_id', 'phenotype_id'])
    if eqtl_df.empty:
        return eqtl_df
    eqtl_df['pval'] = eqtl_df.groupby(  # assign min pval to duplicates
        ['sid', 'pid', 'tissue'])['pval'].transform('min')
    eqtl_df['b'] = eqtl_df.groupby(  # assign min b to duplicaates
        ['sid', 'pid', 'tissue'])['b'].transform('min')
    eqtl_df['b_se'] = eqtl_df.groupby(
        ['sid', 'pid', 'tissue'])['b_se'].transform('min')
    eqtl_df = eqtl_df.drop_duplicates()
    eqtl_df['adj_pval'] = multitest.multipletests(
        eqtl_df['pval'], method='fdr_bh')[1]
    cols = ['sid', 'pid', 'sid_chr', 'sid_pos',
            'adj_pval', 'pval', 'b', 'b_se', 'maf', 'tissue']
    logger.write('  * {} eQTL associations tested'.format(len(eqtl_df)))
    logger.write('  * {} eQTL associations passed FDR <= {}.'.format(
        len(eqtl_df[eqtl_df['adj_pval'] <= fdr_threshold]),
        fdr_threshold))
    logger.write('  * eQTLs mapped at MAF >= {} and pval threshold <={}.'.format(
        maf_threshold, pval_threshold))
    logger.write('  Time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60))
    return eqtl_df[cols]


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        '-i', '--input-files',  nargs='+', default=None,
        help='''Path to the genes.txt file(s) generated by genes.py.''')
    parser.add_argument(
        '-b', '--batches_dir', default=None,
        help='''Filepath to the directory containing batch subdirectories
        'each containing genes.txt.''')
    parser.add_argument(
        "-o", "--output-dir", required=True,
        help="The directory in which to output results.")
    parser.add_argument(
        '--pval-threshold', default=1,
        help='Maximum p value for mapping eQTLs. Default: 1.')
    parser.add_argument(
        '--maf-threshold', default=0.1,
        help='Minimum MAF for variants to include. Default: 0.1.')
    parser.add_argument(
        '-p', '--num-processes', type=int,
        default=int(multiprocessing.cpu_count()/2),
        help='Number of CPUs to use (default: half the number of CPUs).')
    parser.add_argument(
        '-f', '--fdr-threshold', type=float, default=0.05,
        help='''The FDR threshold to consider an eQTL
        statistically significant (default: 0.05).''')
    parser.add_argument(
        '-t', '--tissues', nargs='+',
        help='''Space-separated list of eQTL tissues to query. 
        Note that tissues must be from the same eQTL projects.
        Default is all tissues from the GTEx project.
        User 'codes3d.py --list-eqtl-tissues' for a list of installed tissues.''')
    parser.add_argument(
        '--eqtl-project', type=str, default=None,
        help='''The eQTL project to query. Default: GTEx.
        'use \'codes3d.py --list-eqtl-db\' to list available databases.''')
    parser.add_argument(
        '-c', '--config',
        default=os.path.join(os.path.dirname(
            __file__), '../docs/codes3d.conf'),
        help='The configuration file to use (default: docs/codes3d.conf).')
    parser.add_argument(
        '--suppress-intermediate-files', action='store_true', default=False,
        help='''Do not produce intermediate
        files. These can be used to run the pipeline from 
        an intermediate stage in the event of interruption 
        (default: False).''')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if (not args.input_files) and (not args.batches_dir):
        sys.exit('Required --input-files or --batches_dir missing.')
    c = codes3d.CODES3D(args.config)
    start_time = time.time()
    commons_db = create_engine(c.commons_db_url, echo=False,
                               poolclass=NullPool)
    tissues = codes3d.parse_tissues(
        args.tissues, None, args.eqtl_project, commons_db)
    eqtl_project = tissues['project'].tolist()[0]
    config = configparser.ConfigParser()
    config.read(args.config)
    covariates_dir = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'COVARIATES_DIR'))
    expression_dir = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'EXPRESSION_DIR'))
    eqtl_project_db = create_engine(
        c.eqtl_db_url.format(eqtl_project.lower()), echo=False,
        poolclass=NullPool)
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    if args.input_files:
        gene_files = [os.path.join(os.path.dirname(fp), 'genes.txt')
                      for fp in args.input_files]
    elif args.batches_dir:
        gene_files = [os.path.join(args.batches_dir, batch_dir,  'genes.txt')
                      for batch_dir in os.listdir(args.batches_dir)]
    gene_df = codes3d.parse_intermediate_files(
        gene_files, args.output_dir, 'genes')
    logger = codes3d.Logger(logfile=os.path.join(
        args.output_dir, 'codes3d.log'))
    eqtl_df = map_eqtls(
        gene_df,
        tissues,
        args.num_processes,
        eqtl_project_db,
        covariates_dir,
        expression_dir,
        args.pval_threshold,
        args.maf_threshold,
        args.fdr_threshold,
        logger
    )
    if not args.suppress_intermediate_files:
        gene_df.to_csv(os.path.join(
            args.output_dir, 'genes.txt'), sep='\t', index=False)
    eqtl_df.to_csv(os.path.join(
            args.output_dir, 'eqtls.txt'), sep='\t', index=False)
    logger.write('Done.')
    logger.write('  Total time elasped: {:.2f} mins.'.format(
        (time.time()-start_time)/60))

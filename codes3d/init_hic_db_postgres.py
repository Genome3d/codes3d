#!usr/bin/env python
import os
import pandas as pd
from tqdm import tqdm
import argparse
from sqlalchemy import create_engine
from sqlalchemy_utils import create_database, database_exists
import sys
import codes3d
import time
import gzip
import sqlite3

def create_db(db_auth, db_url, tblspace):
    if tblspace and db_auth.startswith('postgres'):
        db = create_engine(os.path.join(db_auth, 'codes3d'), echo=False)
        sql = '''CREATE DATABASE {} TABLESPACE {}'''
        conn = db.connect()
        conn.execute('commit')
        conn.execute(sql.format(db_url, tblspace))
        conn.close()
    elif db_auth.startswith('sqlite'):
        create_database(os.path.join(db_auth, db_url))
    if database_exists(os.path.join(db_auth, db_url)):
        print('{} database created.'.format(db_url))
    else:
        print('FATAL: {} database not created.'.format(db_url))
        sys.exit()


def build_hic_tables(cell_line, enzyme, fp, db_auth, mapq_cutoff, tblspace):
    start_time = time.time()
    db_url = 'hic_{}_{}'.format(enzyme.lower(), cell_line.lower())
    #db_url = os.path.join(db_auth, 'hic_{}_{}'.format(enzyme.lower(), cell_line.lower()))
    if db_auth.startswith('sqlite'):
        db_url = db_url + '.db'
    if not database_exists(os.path.join(db_auth, db_url)):
        create_db(db_auth, db_url, tblspace)
    db = create_engine(os.path.join(db_auth, db_url), echo=False)
    rep = os.path.basename(fp).split('.')[0]
    rep = rep[:rep.rfind('_merged')].split('_')[-1]
    desc = '  * {} replicate {}'.format(cell_line, rep)
    bar_format = ''
    t = tqdm(total=0, desc=desc, bar_format=bar_format)
    chunksize = 200000
    idx = 1
    #db = create_engine(db_url, echo=False)
    for df in pd.read_csv(
            fp, sep=' ', header=None,
            compression='gzip', chunksize=chunksize):
        '''
        Columns are for short format described 
        https://github.com/aidenlab/juicer/wiki/Pre#medium-format-most-common
        '''
        cols = [ 
            'readname', 'str1', 'chr1', 'pos1', 'frag1', 'str2', 'chr2',
            'pos2', 'frag2', 'mapq1', 'mapq2']
        df.columns = cols
        df = (df.drop(['mapq1', 'mapq2'], axis=1).join( # Remove non-numeric MAPQ as seen in MitDNA
            df[['mapq1', 'mapq2']].apply(pd.to_numeric, errors='coerce')))
        df = df[df[['mapq1', 'mapq2']].notnull().all(axis=1)]
        self_ligation = (  # Filter self-ligating fragments
            (df['chr1'] == df['chr2']) &\
            (df['frag1'] == df['frag2']))
        low_mapq = ( # Filter low quality reads
            (df['mapq1'] >= mapq_cutoff) & (df['mapq2'] >= mapq_cutoff))
        db_cols = ['chr1', 'frag1', 'chr2', 'frag2']
        df = df[~self_ligation & low_mapq][db_cols]
        df_reverse = df.rename( # Reversing to make SQL queries faster
            columns={'chr1': 'chr2',
                     'frag1': 'frag2',
                     'chr2': 'chr1',
                     'frag2': 'frag1'})
        df = pd.concat([
            df.sort_values(by=['chr1', 'frag1']),
            df_reverse.sort_values(by=['chr1', 'frag1'])])
        df[['chr1', 'chr2']] = df[['chr1', 'chr2']].astype(str)
        patched_chrom = (df['chr1'].str.contains('_'))
        df = df[~patched_chrom]
        if_exists = 'replace' if t.total == 0 else 'append'
        df.to_sql(
            rep.lower(), con=db, if_exists=if_exists, index=False)
        t.total += len(df)
        t.update(len(df))
    t.close()
    print('  * Creating index...')
    # '+', '-' throw a sqlite3.OperationalError
    cell_index = cell_line.lower().replace('+', '')
    cell_index = cell_index.replace('-', '_')
    db.execute(
        '''CREATE INDEX idx_{}_{}_chr1_frag1 ON {}
        (chr1, frag1)'''.format(
            cell_index,
            rep.lower(),
            rep.lower()))
    print('    Time elapsed: {:.2f} mins'.format((time.time()-start_time)/60))




def parse_args():
    parser = argparse.ArgumentParser(
        description='Build PostgreSQL database of Hi-C libraries for CoDeS3D.')
    parser.add_argument(
        '--hic-dir', required=True,
        help='''The directory containing (subdirectories containing) 
        Hi-C non_dups.txt.gz files. Directories should be named by the cell lines.''')
    parser.add_argument(
        '-m', '--mapq-cutoff', default=30,
        help='Mininum mapping quality score (mapq) of hic reads. (Default: 30)')
    parser.add_argument(
        "-a", "--db-auth", required=True,
        help='''Connection string of database e.g 'postgresql://codes3d:<password>@localhost'
        or 'sqlite:///<path to write database>' ''')
    parser.add_argument(
        "-e", "--enzyme", required=True,
        help='The restriction enzyme used to prepare the Hi-C library')
    parser.add_argument(
        "-t", "--tablespace", default=None,
        help='Tablespace of database.')
    return parser.parse_args()

    

if __name__ == '__main__':
    args = parse_args()
    start_time = time.time()
    if args.db_auth.startswith('postgres') and not args.tablespace:
        sys.exit('''PostgreSQL databases require a TABLESPACE to be set.
        Provide the name of the tablespace with the '--tablespace' flag.
        See https://www.postgresqltutorial.com/postgresql-create-tablespace/''')
    logfile = '/tmp/codes3d_db.log'
    done = []
    if os.path.exists(logfile):
        f = open(logfile, 'r')
        done = [lib.rstrip() for lib in f]
    print('Building...')    
    for first_level in os.listdir(args.hic_dir):
        if os.path.isdir(os.path.join(args.hic_dir, first_level)):
            for second_level in os.listdir(os.path.join(args.hic_dir, first_level)):
                cell_line = first_level
                fp = os.path.join(args.hic_dir, first_level, second_level)
                if fp in done:
                    continue
                build_hic_tables(cell_line, args.enzyme, fp, args.db_auth, args.mapq_cutoff,
                    args.tablespace)
                with open(logfile, 'a') as f:
                    f.write('{}\n'.format(fp))
        elif os.path.isfile(os.path.join(args.hic_dir, first_level)):
            cell_line = ''
            if args.hic_dir.endswith('/'):
                cell_line = args.hic_dir.strip().split('/')[-2]
            else:
                cell_line = args.hic_dir.strip().split('/')[-1]
            fp = os.path.join(args.hic_dir,first_level)
            if fp in done:
                continue
            build_hic_tables(cell_line, args.enzyme, fp, args.db_auth, args.mapq_cutoff,
                             args.tablespace)
            with open(logfile, 'a') as f:
                f.write('{}\n'.format(fp))

    os.remove(logfile)
    print('Total time elapsed: {:.2f} mins'.format((time.time()-start_time)/60))



#!usr/bin/env python
import os
import pandas as pd
from tqdm import tqdm
import argparse
from sqlalchemy import create_engine
from sqlalchemy_utils import create_database, database_exists
import sys


def build_hic_tables(cell_line, enzyme, fp, db_auth, mapq_cutoff, tblspace):
    db_url = os.path.join(db_auth, 'hic_{}_{}'.format(enzyme.lower(), cell_line.lower()))
    if db_url.startswith('sqlite'):
        db_url = db_url + '.db'
    db = create_engine(db_url, echo=False)
    if not database_exists(db_url):
        create_database(db_url)
        if tblspace and db_url.startswith('postgres'):
            db.execute(
                '''ALTER DATABASE {} SET TABLESPACE {}; '''.format(
                    'hic_{}_{}'.format(db_url), tblspace))
        print('{} database created.'.format(db_url))
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
        idx += len(df)
        t.total += len(df)
        t.update(len(df))
    t.close()
    db.execute(
        '''CREATE INDEX idx_{}_{}_chr1_frag1 ON {}
        (chr1, frag1)'''.format(cell_line, rep.lower(), rep.lower()))



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
                    print("\tProcessed %d%%..." %
                          ((float(i) / float(lines)) * 100))
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



def parse_args():
    parser = argparse.ArgumentParser(
        description='Build SQLite database of  Hi-C libraries for CoDeS3D.')
    parser.add_argument(
        '--hic-dir', required=True,
        help='The directory containing Hi-C non_dups.txt.gz files.')
    parser.add_argument(
        '-m', '--mapq-cutoff', default=30,
        help='Mininum mapping quality score (mapq) of hic reads.')
    parser.add_argument(
        "-a", "--db-auth", required=True,
        help='URL of database e.g postgresql://user:password@hostname')
    parser.add_argument(
        "-e", "--enzyme", required=True,
        help='The restriction enzyme used to prepare the Hi-C library')
    parser.add_argument(
        "-t", "--tablespace", default=None,
        help='Tablespace of database.')
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used for this " +
                        "run (default: codes3d.conf)")

    return parser.parse_args()

    

if __name__ == '__main__':
    args = parse_args()

    print('Building...')
    for first_level in os.listdir(args.hic_dir):
        if os.path.isdir(os.path.join(args.hic_dir, first_level)):
            for second_level in os.listdir(os.path.join(args.hic_dir, first_level)):
                cell_line = first_level
                fp = os.path.join(args.hic_dir,first_level, second_level)
            build_hic_tables(cell_line, args.enzyme, fp, args.db_auth, args.mapq_cutoff, args.tablespace)
        elif os.path.isfile(os.path.join(args.hic_dir, first_level)):
            cell_line = ''
            if args.hic_dir.endswith('/'):
                cell_line = args.hic_dir.strip().split('/')[-2]
            else:
                cell_line = args.hic_dir.strip().split('/')[-1]
            fp = os.path.join(args.hic_dir,first_level)
            build_hic_tables(cell_line, args.enzyme, fp, args.db_auth, args.mapq_cutoff, args.tablespace)
                
                

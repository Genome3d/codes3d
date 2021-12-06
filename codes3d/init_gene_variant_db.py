import pandas as pd
import sys
from tqdm import tqdm
import argparse
import os
import time
from sqlalchemy import create_engine
import codes3d


def build_variants_table(
        variant_fp,
        db_url,
        table):
    """Create a table of genotype. 
    genotype(snp, sample_id..) 
    Indexes: 
    "ix_genotype_snp" btree (snp)
    """
    desc = 'Building variants table...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='variants', desc=desc, disable=False,
             bar_format=bar_format)
    chunksize = 200000
    cols = ['chr', 'frag_id', 'id', 'variant_id']
    db = create_engine(db_url, echo=False)
    idx = 1
    for df in pd.read_csv(
            variant_fp, sep='\t', compression='gzip', chunksize=chunksize):
        df = df.rename(columns={
            'chr': 'chrom',
            'variant_pos': 'locus',
            'rs_id_dbSNP151_GRCh38p7': 'rsid'
            }
        )
        df['id'] = range(idx, idx+len(df))
        cols = ['id', 'rsid', 'chrom', 'locus', 'variant_id', 'ref', 'alt', 'num_alt_per_site']
        if_exists = 'replace' if t.total == 0 else 'append'
        df[cols].to_sql(table, con=db, if_exists=if_exists, index=False)
        idx += len(df)
        t.total += len(df)
        t.update(len(df))
    t.close()
    db.execute('''CREATE INDEX idx_{}_rsid_chrom_locus ON {}
    (rsid, chrom, locus)'''.format(
        '{}_{}'.format(table, 'rsid_chrom_locus'), table))
    db.execute('''CREATE INDEX idx_{}_rsid_chrom_locus ON {}
    (rsid, chrom, locus)'''.format(
        '{}_{}'.format(table, 'id'), table))

    
def build_gene_table(
        gene_fp,
        db_url,
        table):
    """Create a table of gene info. 
    gene(id, chrom, start, end, name, gencode_id) 
    Indexes: 
    "genes_pk" PRIMARY KEY, btree (id)
    "genes_name_id" btree (id)
    "genes_name_index" btree (name)
    """
    desc = 'Building gene table...'
    bar_format = '{desc}: {n_fmt} {unit}'
    t = tqdm(total=0, unit='genes', desc=desc, disable=False,
             bar_format=bar_format)
    chunksize = 2000
    db = create_engine(db_url, echo=False)
    for df in pd.read_csv(
            query_fp, sep='\t', compression='infer', chunksize=chunksize
    ):
        cols = ['frag_id', 'chrom', 'id']
        if_exists = 'replace' if t.total == 0 else 'append'
        df[cols].to_sql(table, con=db, if_exists=if_exists, index=False)
        t.total = len(df)
        t.update(len(df))
    t.close()
    create_index(table, db)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Build database tables gene and variant fragments')
    parser.add_argument(
        '-v', '--variants',
        help='The filepath to the variant lookup file')
    parser.add_argument(
        '-g', '--genes',
        help='''The filepath to a sorted gene reference file containing
              'id, chrom, start, end, name, and gencode_id columns'.
              Note that the id indicates the sorting (chrom, start)
              order of the file starting at 1.''')
    parser.add_argument(
        '-t', '--table',
        help='Name of table e.g. variant_lookup_mboi.')
    parser.add_argument(
        '-u', '--db-url', required=True,
        help='''URL of database e.g postgresql://user:password@hostname/database
        This is usually codes3d_commons''')
    parser.add_argument(
        '-c', '--config',
        default=os.path.join(os.path.dirname(__file__),
                             '../docs/codes3d.conf'),
        help='''The configuration file to be used for this
        run (default: codes3d.conf)''')

    args = parser.parse_args()
    c = codes3d.CODES3D(args.config)
    start_time = time.time()
    if not args.variants and not args.genes:
        sys.exit('Variants or gene bed file missing.')
    if args.variants:
        build_variants_table(args.variants, args.db_url, args.table)
    elif args.genes:
        build_gene_table(args.genes, args.db_url, args.table)
    print('Time elapsed: {:.2f} mins'.format(
        (time.time() - start_time)/60))

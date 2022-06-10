#! /usr/bin/env python

import pandas as pd
import numpy as np
import multiprocessing
import tqdm
import sys
import os
import sqlite3
from sqlalchemy import create_engine
import time
from itertools import repeat


def get_cell_interactions(
        spatial,
        enzyme,
        replicate,
        query_type,
        df,
        interactions):
    cell_line = os.path.dirname(replicate).split('/')[-1]
    rep_interactions = []
    
    # TODO: Add implementation for PostgreSQL database
    
    if spatial == 'hic':
        sql = '''
              SELECT chr2, fragment2 FROM interactions WHERE chr1='{}' AND fragment1={}
              '''
        with create_engine('sqlite:///{}'.format(replicate), echo=False).connect() as con:
            for idx, row in df.iterrows():
                from_db = pd.read_sql_query(sql.format(
                    row['chrom'][3:], row['fragment']), con=con)
                if not from_db.empty:
                    from_db = from_db.rename(
                        columns={'chr2': 'fragment_chr', 'fragment2': 'fragment'})
                    #from_db[query_type] = row[query_type]
                    from_db['fragment_chr'] = 'chr' + \
                        from_db['fragment_chr'].astype(str)
                    from_db['query_type'] = query_type
                    from_db['query_chr'] = row['chrom']
                    from_db['query_fragment'] = row['fragment']
                    from_db['fragment'] = from_db['fragment'].astype('int64')
                    rep_interactions.append(from_db)
        if len(rep_interactions) == 0:
            return
        rep_interactions = pd.concat(rep_interactions)
        # A fragment can have many interactions with the same fragment
        rep_interactions = rep_interactions.drop_duplicates()  # Keep unique replicates
        rep = os.path.basename(replicate).split('.')[0]
        rep = rep[:rep.rfind('_merged')].split('_')[-1]
        rep_interactions['replicate'] = rep
        rep_interactions['cell_line'] = cell_line
        rep_interactions['enzyme'] = enzyme
        self_ligation = (  # Filter self-ligating fragments
            (rep_interactions['fragment_chr'] == rep_interactions['query_chr']) &\
            (rep_interactions['fragment'] == rep_interactions['query_fragment']))
        rep_interactions = rep_interactions[~self_ligation]
        interactions.append(rep_interactions)

    else:
        sql = '''SELECT * from interactions WHERE p_fid={} OR oe_fid={}'''
        with create_engine('sqlite:///{}'.format(replicate), echo=False).connect() as con:
            for idx, row in df.iterrows():
                from_db = pd.read_sql_query(sql.format(
                    row['fragment'], row['fragment']), con=con)
                if not from_db.empty:
                    from_db_subset = from_db[['p_fid','p_name','oe_fid','oe_name',
                        'n_reads','score']]
                    from_db_subset['query_type'] = query_type
                    from_db_subset['query_fragment'] = row['fragment']
                    rep_interactions.append(from_db_subset)
        if len(rep_interactions) == 0:
            return
        rep_interactions = pd.concat(rep_interactions)
        rep_interactions = rep_interactions.drop_duplicates()
        rep = os.path.basename(replicate).split('.')[0]
        rep_interactions['replicate'] = rep
        rep_interactions['cell_line'] = cell_line
        rep_interactions['enzyme'] = enzyme
        interactions.append(rep_interactions)

def filter_inter_for_pchic(interactions):
    
    pchic_interactions = []
    for row in interactions.itertuples(index=False):
        # remove snps that interact with unannotated oe_frag
        if row.query_fragment == row.p_fid  and \
                row.oe_name != ".":
            #row['snp_in_promoter'] = row['query_snp']
            pchic_interactions.append(row)
        elif row.query_fragment == row.oe_fid and \
                row.oe_name != ".":
            pchic_interactions.append(row)
        elif row.query_fragment == row.oe_fid and \
                row.oe_name == ".":
            pchic_interactions.append(row)
        elif row.query_fragment == row.p_fid and \
                row.oe_name == ".":
            pass
    pchic_interactions = pd.DataFrame(pchic_interactions)
    
    return pchic_interactions

def find_interactions(
        spatial,
        fragment_df,
        query_type,
        lib_dir,
        hic_df,
        pchic_df,
        num_processes,
        logger,
        suppress_intermediate_files=False):
    """Finds fragments that interact with gene fragments.
    """
    start_time = time.time()
    logger.write("Finding chromatin interactions in...")
    restriction_enzymes = fragment_df['enzyme'].drop_duplicates().tolist()
    interactions_df = []
    for enzyme in restriction_enzymes:
        if spatial == 'hic':
            _3dgi_data_dir = os.path.join(lib_dir, 'hic', enzyme, 'hic_data')
        else:
            _3dgi_data_dir = os.path.join(lib_dir, 'pchic', enzyme, 'pchic_data')
        manager = multiprocessing.Manager()
        gene_interactions = manager.list()
        enzyme_qc_df = manager.list()
        enzyme_qc_all_df = manager.list()
        num_processes = int(min(16, multiprocessing.cpu_count()/2))
        if spatial == 'hic':
            cell_lines = hic_df[hic_df['enzyme'] == enzyme]['library'].tolist()
        else:
            cell_lines = pchic_df[pchic_df['enzyme'] == enzyme]['library'].tolist()
        replicates = []
        for cell_line in cell_lines:  # TODO: Replace with database
            replicates += [os.path.join(_3dgi_data_dir, cell_line, replicate)
                           for replicate in os.listdir(
                                            os.path.join(_3dgi_data_dir, cell_line))
                                            if replicate.endswith('.db')]
        enzyme_df = fragment_df[fragment_df['enzyme'] == enzyme]
        enzyme_df = enzyme_df[['fragment',
                               'chrom', 'enzyme']].drop_duplicates()
        with multiprocessing.Pool(num_processes) as pool:
            if spatial == 'pchic':
                desc = '  * PCHi-C libraries restricted with {}'.format(enzyme)
            else:
                desc = '  * Hi-C libraries restricted with {}'.format(enzyme)
            bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
            bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
            for _ in tqdm.tqdm(
                    pool.istarmap(get_cell_interactions,
                                  zip(
                                      repeat(spatial),
                                      repeat(enzyme),
                                      replicates,
                                      repeat(query_type),
                                      repeat(enzyme_df),
                                      repeat(gene_interactions)
                                  )
                                  ),
                    total=len(replicates),
                    desc=desc,
                    unit='batches',
                    ncols=80,
                    bar_format=bar_format):
                pass
        for df in gene_interactions:
            interactions_df.append(df)
    if not interactions_df:
        msg = '''Program exiting: Chromatin interactions involving input
                SNP(s)/gene(s) could not be found.'''
        sys.exit(msg)
        logger.write(msg)
    else:
        interactions_df = pd.concat(interactions_df)
    
    if spatial == 'hic':
        cols = ['cell_line', 'query_type', 'query_chr', 'query_fragment',
                'fragment_chr', 'fragment', 'replicate', 'enzyme']
        interactions_df = interactions_df[cols]
        interactions_df = interactions_df.drop_duplicates()
        logger.write('  Time elasped: {:.2f} mins.'.format((time.time() - start_time)/60))

    else:
        cols = ['p_fid', 'p_name', 'oe_fid', 'oe_name', 'n_reads',
                'score', 'query_type', 'query_fragment', 
                'replicate', 'cell_line', 'enzyme']
        interactions_df = interactions_df[cols]
        interactions_df = interactions_df.drop_duplicates()
        interactions_df = filter_inter_for_pchic(interactions_df)
        logger.write('  Time elasped: {:.2f} mins.'.format((time.time() - start_time)/60))

    return interactions_df

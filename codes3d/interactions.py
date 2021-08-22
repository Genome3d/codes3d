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
        enzyme,
        replicate,
        query_type,
        df,
        interactions,
        lib_dir,
):
    cell_line = os.path.dirname(replicate).split('/')[-1]
    rep_interactions = []
    
    # TODO: Add implementation for PostgreSQL database
    
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


def find_interactions(
        fragment_df,
        query_type,
        lib_dir,
        hic_df,
        #output_dir,
        num_processes,
        logger,
        suppress_intermediate_files=False):
    """Finds fragments that interact with gene fragments.
    """
    start_time = time.time()
    logger.write("Finding chromatin interactions in...")
    # print(fragment_df)
    restriction_enzymes = fragment_df['enzyme'].drop_duplicates().tolist()
    interactions_df = []
    for enzyme in restriction_enzymes:
        hic_data_dir = os.path.join(lib_dir, 'hic', enzyme, 'hic_data')
        manager = multiprocessing.Manager()
        gene_interactions = manager.list()
        enzyme_qc_df = manager.list()
        enzyme_qc_all_df = manager.list()
        num_processes = int(min(16, multiprocessing.cpu_count()/2))
        cell_lines = hic_df[hic_df['enzyme'] == enzyme]['library'].tolist()
        replicates = []
        for cell_line in cell_lines:  # TODO: Replace with database
            replicates += [os.path.join(hic_data_dir, cell_line, replicate)
                           for replicate in os.listdir(
                os.path.join(hic_data_dir, cell_line))
                if replicate.endswith('.db')]
        enzyme_df = fragment_df[fragment_df['enzyme'] == enzyme]
        enzyme_df = enzyme_df[['fragment',
                               'chrom', 'enzyme']].drop_duplicates()

        with multiprocessing.Pool(num_processes) as pool:
            desc = '  * Hi-C libraries restricted with {}'.format(enzyme)
            bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
            bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
            for _ in tqdm.tqdm(
                    pool.istarmap(get_cell_interactions,
                                  zip(
                                      repeat(enzyme),
                                      replicates,
                                      repeat(query_type),
                                      repeat(enzyme_df),
                                      repeat(gene_interactions),
                                      repeat(lib_dir),
                                      #repeat(enzyme_qc_df),
                                      #repeat(enzyme_qc_all_df)
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
    interactions_df = pd.concat(interactions_df)
    cols = ['cell_line', 'query_type', 'query_chr', 'query_fragment',
            'fragment_chr', 'fragment',  # 'fragment_distance',
            'replicate', 'enzyme']
    interactions_df = interactions_df[cols]
    interactions_df = interactions_df.drop_duplicates()
    logger.write('  Time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60))
    return interactions_df

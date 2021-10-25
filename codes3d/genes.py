#! /usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import time
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import multiprocessing
import tqdm
from itertools import repeat


def get_gene_fragments(gene_df, restriction_enzymes, db):
    db.dispose()
    gene_df = gene_df.sort_values(by=['id'])
    fragment_df = []
    chunksize = 1000
    chunks = [gene_df[i:i+chunksize]
              for i in range(0, gene_df.shape[0], chunksize)]
    for enzyme in restriction_enzymes:
        table = 'gene_lookup_{}'
        if enzyme in ['MboI', 'DpnII']:  # MboI and DpnII have the same restriction sites
            table = table.format('mboi')
        else:
            table = table.format(enzyme.lower())
        sql = '''SELECT * FROM {} WHERE id >= {} AND id <= {}'''
        with db.connect() as con:
            for chunk in chunks:
                df = pd.read_sql(sql.format(
                    table, chunk['id'].min(), chunk['id'].max()), con=con)
                df['enzyme'] = enzyme
                fragment_df.append(df)
    fragment_df = pd.concat(fragment_df)
    gene_df = pd.merge(gene_df, fragment_df, how='inner', on=['id', 'chrom'])
    return gene_df


def process_snp_genes(
        gene_info_df,
        inter_df,
        snp_df,
        gene_fp):
    inter_df = inter_df.rename(columns={'fragment_chr': 'chr'})
    gene_cols = [
        'snp', 'chr', 'locus', 'variant_id',
        'gene', 'gencode_id', 'gene_chr', 'gene_start', 'gene_end',
        'interactions', 'replicates', 'cell_line', 'cell_line_replicates',
        'num_cell_lines']
    gene_info_df = gene_info_df.rename(columns={
        'name': 'gene',
        'chr': 'gene_chr',
        'start': 'gene_start',
        'end': 'gene_end'})
    gene_df = snp_df.merge(
        inter_df, how='inner',
        on=['chr', 'fragment', 'enzyme'])
    gene_df = gene_df.merge(
        gene_info_df, how='left',
        left_on=['query_chr', 'query_fragment', 'enzyme'],
        right_on=['gene_chr', 'gene_fragment', 'enzyme'])
    gene_df = gene_df[gene_cols].drop_duplicates()
    gene_df.to_csv(
        gene_fp, sep='\t', mode='a', header=False, index=False)
    return gene_df



def find_snp_genes(
        chunk_df,
        enzyme,
        enzyme_genes):
    db.dispose()
    chunk_df = chunk_df.sort_values(by=['fragment'])
    chrom = chunk_df['fragment_chr'].drop_duplicates().to_list()[0]
    table = 'gene_lookup_{}'
    if enzyme in ['MboI', 'DpnII']:  # MboI and DpnII have the same restriction sites
        table = table.format('mboi')
    else:
        table = table.format(enzyme.lower())
    chunk_df['fragment'] = chunk_df['fragment'].astype(int)
    sql = ''' SELECT * FROM {0}
    JOIN genes on {0}.id=genes.id
    WHERE {0}.chrom = '{1}' AND {0}.frag_id >= {2} AND {0}.frag_id <= {3}'''
    sql = sql.format(table, chrom,
                     chunk_df['fragment'].min(), chunk_df['fragment'].max())
    df = pd.DataFrame()
    with db.connect() as con:
        df = pd.read_sql_query(sql, con)
    if df.empty:
        return
    df = df.loc[:, ~df.columns.duplicated()]
    df = df.rename(
        columns={'id': 'gene_id', 'name': 'gene', 'chrom': 'gene_chr',
                 'start': 'gene_start', 'end': 'gene_end'})
    chunk_df = pd.merge(chunk_df, df, how='inner',
                        left_on=['fragment_chr', 'fragment'], right_on=['gene_chr', 'frag_id'])

    enzyme_genes.append(chunk_df)


def fetch_hic_libs(db):
    with db.connect() as con:
        hic_libs = pd.read_sql_query(
            'SELECT library, enzyme, rep_count FROM meta_hic',
        con=con)
        return hic_libs.drop_duplicates()


def get_gene_by_id(
        snp_df,
        inter_df,
        _db,
        logger,
):
    logger.write('Identifying genes interacting with SNPs in...')
    global db
    db = _db
    start_time = time.time()
    enzymes = inter_df['enzyme'].drop_duplicates().tolist()
    all_genes_df = []
    #db = create_engine(db_url, echo=False, poolclass=NullPool)
    hic_libs = fetch_hic_libs(db)
    hic_libs = hic_libs.rename(columns={'rep_count': 'cell_line_replicates'})
    for enzyme in enzymes:
        manager = multiprocessing.Manager()
        num_processes = int(min(16, multiprocessing.cpu_count()/2))
        enzyme_genes = manager.list()
        enzyme_df = []
        with multiprocessing.Pool(processes=num_processes) as pool:
            df = inter_df[inter_df['enzyme'] == enzyme]
            snp_interactions = [
                df[df['fragment_chr'] == chrom]
                for chrom in df['fragment_chr'].drop_duplicates().to_list()
            ]
            desc = '  * Hi-C libraries restricted with {}'.format(enzyme)
            bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
            '''
            for snp in snp_interactions:
                find_snp_genes(snp,
                               enzyme,
                               enzyme_genes,
                               db)
            '''
            for _ in tqdm.tqdm(pool.istarmap(
                    find_snp_genes,
                    zip(snp_interactions,
                        repeat(enzyme),
                        repeat(enzyme_genes))),
                    total=len(snp_interactions), desc=desc, unit='batches',
                    ncols=80, bar_format=bar_format):
                pass
            
        for df in enzyme_genes:
            enzyme_df.append(df)
        enzyme_df = pd.concat(enzyme_df)
        enzyme_df = enzyme_df.merge(
            hic_libs, how='left',
            left_on=['cell_line', 'enzyme'], right_on=['library', 'enzyme'])
        enzyme_df['interactions'] = enzyme_df.groupby(
            ['query_chr', 'query_fragment', 'gencode_id', 'cell_line'])[
            'fragment'].transform('count')
        df = enzyme_df[
            ['query_chr', 'query_fragment', 'gencode_id', 'cell_line', 'replicate']
        ].drop_duplicates()
        df['replicates'] = df.groupby(
            ['query_chr', 'query_fragment', 'gencode_id', 'cell_line'])[
            'replicate'].transform('count')
        enzyme_df = enzyme_df.merge(
            df, how='left',
            on=['query_chr', 'query_fragment', 'gencode_id', 'cell_line', 'replicate'])
        enzyme_df = enzyme_df.drop(
            columns=['fragment_chr', 'fragment',
                     'frag_id', 'library', 'replicate']
        )
        enzyme_df = enzyme_df.drop_duplicates()
        all_genes_df.append(enzyme_df.drop_duplicates())
    logger.write('  * Filtering gene interactions...')
    all_genes_df = pd.concat(all_genes_df)
    all_genes_df = all_genes_df.merge(
        snp_df, left_on=['query_fragment', 'query_chr', 'enzyme'],
        right_on=['fragment', 'chrom', 'enzyme'],  how='inner'
    )  # .drop_duplicates()
    all_genes_df['sum_cell_lines'] = all_genes_df.groupby(
        ['snp', 'gene'])['cell_line'].transform('count')
    all_genes_df['sum_interactions'] = all_genes_df.groupby(
        ['snp', 'gene'])[
        'interactions'].transform('sum')
    all_genes_df['sum_replicates'] = all_genes_df.groupby(
        ['snp', 'gene'])[
        'replicates'].transform('sum')
    condition = (
        (all_genes_df['interactions'] / all_genes_df['cell_line_replicates'] <= 1) &
        (all_genes_df['sum_replicates'] < 2) &
        (all_genes_df['sum_cell_lines'] < 2))
    gene_df = all_genes_df[~condition].drop_duplicates()
    cols = ['snp', 'chrom', 'locus', 'variant_id',
            'gene', 'gencode_id', 'gene_chr', 'gene_start', 'gene_end',
            'interactions', 'replicates', 'cell_line', 'cell_line_replicates',
            'sum_interactions', 'sum_replicates', 'sum_cell_lines']
    logger.write('  Time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60))
    return gene_df[cols].drop_duplicates()


def get_gene_by_position(df, db):
    gene_df = []
    omitted_genes = []
    sql = '''SELECT * FROM genes where chrom = '{}' and start= {} and "end" = {};'''
    with db.connect() as con:
        for idx, row in df.iterrows():
            try:
                chrom = row['gene'].split(':')[0]
                positions = row['gene'].split(':')[1].split('-')
                res = pd.read_sql_query(
                    sql.format(chrom, positions[0], positions[1]),
                    con=con)
                if res.empty:
                    omitted_genes.append(row['gene'])
                else:
                    gene_df.append(res)
            except IndexError:
                sys.exit(
                    'Exiting: {} not in <chr:start-end> format'.format(row['gene']))
    return gene_df, omitted_genes


def get_gene_by_gencode(df, db):
    gene_df = []
    omitted_genes = []
    sql = '''SELECT * FROM genes where gencode_id LIKE '{}.%%';'''
    with db.connect() as con:
        for idx, row in df.iterrows():
            gencode_id = row['gene'].split('.')[0]
            res = pd.read_sql_query(sql.format(gencode_id), con=con)
            if res.empty:
                omitted_genes.append(row['gene'])
            else:
                gene_df.append(res)
    return gene_df, omitted_genes


def get_gene_by_symbol(df, db):
    gene_df = []
    omitted_genes = []
    sql = '''SELECT * FROM genes where name = '{}';'''
    with db.connect() as con:
        for idx, row in df.iterrows():
            symbol = row['gene'] if 'orf' in row['gene']  else row['gene'].upper() #ORFs are not capitalized.
            res = pd.read_sql_query(sql.format(symbol), con=db)
            if res.empty:
                omitted_genes.append(row['gene'])
            else:
                gene_df.append(res)
    return gene_df, omitted_genes


def get_file_gene_info(df, db):
    gene_df = []
    omitted_genes = []
    df = df[~df['gene'].isna()]
    position_df = df[(df['gene'].str.startswith('chr')) &
                     (df['gene'].str.contains(':'))]
    gencode_df = df[(df['gene'].str.startswith('ENSG')) &
                    (df['gene'].str.len() >= 15)]
    symbol_df = df[~df.isin(position_df)].dropna()
    symbol_df = symbol_df[~symbol_df.isin(gencode_df)].dropna()
    if not position_df.empty:
        temp_gene_df, temp_omitted_genes = get_gene_by_position(
            position_df, db)
        gene_df += temp_gene_df
        omitted_genes += temp_omitted_genes
    if not gencode_df.empty:
        temp_gene_df, temp_omitted_genes = get_gene_by_gencode(gencode_df, db)
        gene_df += temp_gene_df
        omitted_genes += temp_omitted_genes
    if not symbol_df.empty:
        temp_gene_df, temp_omitted_genes = get_gene_by_symbol(symbol_df, db)
        gene_df += temp_gene_df
        omitted_genes += temp_omitted_genes
    return gene_df, omitted_genes


def get_gene_info(
        inputs,
        hic_df,
        output_dir,
        db,
        logger,
        suppress_intermediate_files=False):
    enzymes = hic_df['enzyme'].drop_duplicates().tolist()
    gene_df = []
    omitted_genes = []
    gene_inputs = []
    for inp in inputs:
        inp = inp.strip()
        if os.path.isfile(inp):
            df = []
            with open(inp) as f:
                for line in f:
                    if not line.strip() == '':
                        df.append([line.strip()])
            df = pd.DataFrame(df, columns=['gene'])
            if len(df) > 1:
                #sys.exit(
                #    'Exiting: CoDeS3D currently supports only one gene input')
                pass
            temp_gene_df, temp_omitted_genes = get_file_gene_info(df, db)
            gene_df += temp_gene_df
            omitted_genes += temp_omitted_genes
        else:
            if len(gene_inputs) >= 1:
                #sys.exit(
                #    'Exiting: CoDeS3D currently supports only one gene input')
                pass
            gene_inputs.append(inp)
            if inp.startswith('chr') and ':' in inp:
                df = pd.DataFrame([inp], columns=['gene'])
                temp_gene_df, temp_omitted_genes = get_gene_by_position(df, db)
                gene_df += temp_gene_df
                omitted_genes += temp_omitted_genes
            elif inp.startswith('ENSG') and len(inp) >= 15:
                df = pd.DataFrame([inp], columns=['gene'])
                temp_gene_df, temp_omitted_genes = get_gene_by_gencode(df, db)
                gene_df += temp_gene_df
                omitted_genes += temp_omitted_genes
            else:
                df = pd.DataFrame([inp], columns=['gene'])
                temp_gene_df, temp_omitted_genes = get_gene_by_symbol(df, db)
                gene_df += temp_gene_df
                omitted_genes += temp_omitted_genes

    if len(gene_df) == 0:
        sys.exit('Exiting: We could not find your genes in our databases.')
    else:
        gene_df = pd.concat(gene_df).drop_duplicates()
    if len(omitted_genes) > 0:
        msg = '  * Warning: {} gene(s) not found in our database.'.format(
            len(omitted_genes))
        msg = msg + ':\n\t{}'.format(', '.join(omitted_genes))
        logger.write(msg)

    gene_df = get_gene_fragments(gene_df, enzymes, db)
    gene_df = gene_df.rename(columns={'frag_id': 'fragment'})
    return(gene_df)

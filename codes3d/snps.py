#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
import time
from sqlalchemy import create_engine
from itertools import repeat
import multiprocessing
import tqdm
import genes
import eqtls


def find_snps_hic(
        spatial,
        inter_df,
        gene_info_df,
        tissues,
        output_dir,
        C,
        genotypes_fp,
        _eqtl_project_db,
        covariates_dir,
        expression_dir,
        pval_threshold,
        maf_threshold,
        fdr_threshold,
        num_processes,
        _db,
        logger,
        suppress_intermediate_files=False
):
    start_time = time.time()
    global eqtl_project_db
    eqtl_project_db = _eqtl_project_db
    global db
    db = _db
    enzymes = inter_df['enzyme'].drop_duplicates().tolist()
    _3dgi_libs = genes.fetch_3dgi_libs(spatial, db)
    _3dgi_libs = _3dgi_libs.rename(columns={'rep_count': 'cell_line_replicates'})
    inter_df = inter_df.merge(
            _3dgi_libs, how='left',
            left_on=['cell_line', 'enzyme'], right_on=['library', 'enzyme'])
    default_chrom = ['chr' + str(i)
                     for i in list(range(1, 23))] + ['X', 'Y', 'M']
    chrom_list = inter_df['fragment_chr'].drop_duplicates().tolist()
    chrom_list = [i for i in default_chrom if i in chrom_list]
    inter_df = inter_df[inter_df['fragment_chr'].isin(default_chrom)]
    inter_df = inter_df.astype({'fragment': int})
    gene_info_df = gene_info_df.rename(
        columns={
            'name': 'gene',
            'chr': 'gene_chr',
            'start': 'gene_start',
            'end': 'gene_end',
            'fragment': 'gene_fragment',
            'id': 'gene_id'})
    all_snps = []
    all_genes = []
    all_eqtls = []
    logger.write('Finding SNPs within fragments interacting with genes in...')
    for chrom in sorted(chrom_list):
        chrom_dir = os.path.join(output_dir, chrom)
        #if os.path.exists(os.path.join(chrom_dir, 'eqtls.txt')):
        #    logger.write('  Warning: {} already exists. Skipping.'.format(
        #        os.path.join(chrom_dir, 'eqtls.txt')))
        #    continue
        logger.write(' Chromosome {}'.format(chrom))
        snp_cols = ['snp', 'variant_id', 'chr',
                    'locus', 'id', 'fragment', 'enzyme']
        chrom_df = inter_df[inter_df['fragment_chr'] == chrom]
        chrom_df = chrom_df.astype({'fragment': int,
                                    'fragment_chr': object})
        enzymes = chrom_df['enzyme'].drop_duplicates().tolist()
        snp_df = []
        for enzyme in enzymes:
            enzyme_df = chrom_df[chrom_df['enzyme'] == enzyme]
            enzyme_df = enzyme_df.merge(
                gene_info_df, how='inner',
                left_on=['query_chr', 'query_fragment', 'enzyme'],
                right_on=['chrom', 'gene_fragment', 'enzyme'])
            fragment_df = enzyme_df[
                ['gencode_id', 'fragment_chr', 'fragment']].drop_duplicates()
            enzyme_df = enzyme_df.sort_values(by=['fragment'])
            chunksize = 20000
            enzyme_chunks = [enzyme_df[i:i+chunksize]
                             for i in range(0, enzyme_df.shape[0], chunksize)]
            manager = multiprocessing.Manager()
            snps = manager.list()
            desc = '  * Hi-C libraries restricted with {}'.format(
                enzyme)
            bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
            '''
            for df in tqdm.tqdm(enzyme_chunks, desc=desc, unit='batches',
                                ncols=80, bar_format=bar_format):
                find_gene_snps(
                    df,
                    enzyme,
                    snps)
            '''
            with multiprocessing.Pool(processes=4) as pool:
                for _ in tqdm.tqdm(
                        pool.istarmap(
                            find_gene_snps,
                            zip(enzyme_chunks,
                                repeat(spatial),
                                repeat(enzyme),
                                repeat(snps))
                        ),
                        total=len(enzyme_chunks), desc=desc, unit='batches',
                        ncols=80, bar_format=bar_format):
                    pass

            for df in snps:
                df['enzyme'] = enzyme
                snp_df.append(df)
        if len(snp_df) == 0:
            continue
        snp_df = pd.concat(snp_df)
        logger.verbose = False
        gene_df, snp_df = filter_snp_fragments_hic(
            snp_df, logger)
        snp_df.sort_values(by=['variant_id'], inplace=True)
        snp_list = snp_df['variant_id'].drop_duplicates().tolist()
        batchsize = 2000
        snp_batches = [snp_list[i:i + batchsize]
                       for i in range(0, len(snp_list), batchsize)]
        chrom_eqtl_df = []
        for batch_num, snp_batch in enumerate(snp_batches):
            if len(snp_batches) > 1:
                logger.verbose = True
                logger.write('  Mapping eQTLs batch {} of {}'.format(
                    batch_num+1, len(snp_batches)))
                logger.verbose = False
            batch_gene_df = gene_df[gene_df['variant_id'].isin(snp_batch)]
            eqtl_df = eqtls.map_eqtls(
                batch_gene_df,
                tissues,
                output_dir,
                C,
                genotypes_fp,
                num_processes,
                eqtl_project_db,
                covariates_dir,
                expression_dir,
                pval_threshold,
                maf_threshold,
                fdr_threshold,
                logger)
            if eqtl_df is None:
                continue
            chrom_eqtl_df.append(eqtl_df)
        if len(chrom_eqtl_df) > 0:
            chrom_eqtl_df = pd.concat(chrom_eqtl_df)
        else:
            chrom_eqtl_df = pd.DataFrame()
        if not suppress_intermediate_files:
            os.makedirs(chrom_dir, exist_ok=True)
            snp_df.to_csv(os.path.join(chrom_dir, 'snps.txt'),
                          sep='\t', index=False)
            gene_df.to_csv(os.path.join(chrom_dir, 'genes.txt'),
                           sep='\t', index=False)
            chrom_eqtl_df.to_csv(os.path.join(chrom_dir, 'eqtls.txt'),
                                 sep='\t', index=False)
        all_eqtls.append(chrom_eqtl_df)
        all_snps.append(snp_df)
        all_genes.append(gene_df)
        logger.verbose = True
    if len(all_eqtls) == 0:
        snp_df = pd.DataFrame()
        gene_df = pd.DataFrame()
        eqtl_df = pd.DataFrame()
    else:
        snp_df = pd.concat(all_snps)
        gene_df = pd.concat(all_genes)
        eqtl_df = pd.concat(all_eqtls)
    logger.verbose = True
    logger.write('  Time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60))
    return snp_df, gene_df, eqtl_df

def find_snps_pchic(
        spatial,
        inter_df,
        gene_info_df,
        tissues,
        output_dir,
        C,
        genotypes_fp,
        _eqtl_project_db,
        covariates_dir,
        expression_dir,
        pval_threshold,
        maf_threshold,
        fdr_threshold,
        num_processes,
        _db,
        logger,
        suppress_intermediate_files=False):
    start_time = time.time()
    global eqtl_project_db
    eqtl_project_db = _eqtl_project_db
    global db
    db = _db
    enzymes = inter_df['enzyme'].drop_duplicates().tolist()
    _3dgi_libs = genes.fetch_3dgi_libs(spatial, db)
    _3dgi_libs = _3dgi_libs.rename(columns={'rep_count': 'cell_line_replicates'})
    inter_df = inter_df.merge(
            _3dgi_libs, how='left',
            left_on=['cell_line', 'enzyme'], right_on=['library', 'enzyme'])
    inter_df = inter_df[['p_fid','oe_fid','n_reads','score',
        'query_type', 'query_fragment', 'replicate', 'cell_line', 'enzyme',
        'library']]
    inter_df['inter_frag'] = np.where(inter_df['query_fragment'] == 
                                inter_df['p_fid'], inter_df['oe_fid'], inter_df['p_fid'])
    inter_df = inter_df[['n_reads','score','query_type','query_fragment','inter_frag',
        'replicate','cell_line','enzyme']].drop_duplicates()
    gene_info_df = gene_info_df.rename(
            columns={
                'name': 'gene',
                'chr': 'gene_chr',
                'start': 'gene_start',
                'end': 'gene_end',
                'fragment': 'gene_fragment',
                'id': 'gene_id'})
    all_snps = []
    all_genes = []
    all_eqtls = []
    logger.write('Finding SNPs within fragments interacting with gene promoters in...')
    snp_df = []
    for enzyme in enzymes:
        enzyme_dir = os.path.join(output_dir, enzyme)
        enzyme_df = inter_df[inter_df['enzyme'] == enzyme]
        enzyme_df = enzyme_df.merge(
                gene_info_df, how='inner',
                left_on = ['query_fragment', 'enzyme'],
                right_on = ['gene_fragment', 'enzyme'])
        fragment_df = enzyme_df[
                ['gencode_id','query_fragment','inter_frag','project']].drop_duplicates()
        enzyme_df = enzyme_df.drop(columns=['project']).drop_duplicates()
        enzyme_df = enzyme_df.sort_values(by=['inter_frag'])
        chunksize = 10000
        enzyme_chunks = [enzyme_df[i:i+chunksize]
                for i in range(0, enzyme_df.shape[0], chunksize)]
        manager = multiprocessing.Manager()
        snps = manager.list()
        desc = '  * PCHi-C libraries restricted with {}'.format(
            enzyme)
        bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
        with multiprocessing.Pool(processes=4) as pool:
            for _ in tqdm.tqdm(
                    pool.istarmap(
                        find_gene_snps,
                        zip(enzyme_chunks,
                            repeat(spatial),
                            repeat(enzyme),
                            repeat(snps))
                        ),
                    total=len(enzyme_chunks), desc=desc, unit='batches',
                    ncols=80, bar_format=bar_format): 
                pass
        for df in snps:
            df['enzyme'] = enzyme
            snp_df.append(df)
        if len(snp_df) == 0:
            continue
    snp_df = pd.concat(snp_df)
    logger.verbose = False
    gene_df, snp_df = filter_snp_fragments_pchic(
            snp_df, logger)
    snp_df.sort_values(by=['variant_id'], inplace=True)
    snp_list = snp_df['variant_id'].drop_duplicates().tolist()
    batchsize = 2000
    snp_batches = [snp_list[i:i + batchsize]
            for i in range(0, len(snp_list), batchsize)]
    chrom_eqtl_df = []
    for batch_num, snp_batch in enumerate(snp_batches):
        if len(snp_batches) > 1:
            logger.verbose = True
            logger.write('  Mapping eQTLs batch {} of {}'.format(
                batch_num+1, len(snp_batches)))
            logger.verbose = False
        batch_gene_df = gene_df[gene_df['variant_id'].isin(snp_batch)]
        eqtl_df = eqtls.map_eqtls(
                batch_gene_df,
                tissues,
                output_dir,
                C,
                genotypes_fp,
                num_processes,
                eqtl_project_db,
                covariates_dir,
                expression_dir,
                pval_threshold,
                maf_threshold,
                fdr_threshold,
                logger)
        if eqtl_df is None:
            continue
        chrom_eqtl_df.append(eqtl_df)
    if len(chrom_eqtl_df) > 0:
        chrom_eqtl_df = pd.concat(chrom_eqtl_df)
    else:
        chrom_eqtl_df = pd.DataFrame()
    if not suppress_intermediate_files:
        os.makedirs(enzyme_dir, exist_ok=True)
        snp_df.to_csv(os.path.join(enzyme_dir, 'snps.txt'),
                sep='\t', index=False)
        gene_df.to_csv(os.path.join(enzyme_dir, 'genes.txt'),
                sep='\t', index=False)
        chrom_eqtl_df.to_csv(os.path.join(enzyme_dir, 'eqtls.txt'),
                sep='\t', index=False)
    all_eqtls.append(chrom_eqtl_df)
    all_snps.append(snp_df)
    all_genes.append(gene_df)
    logger.verbose = True
    if len(all_eqtls) == 0:
        snp_df = pd.DataFrame()
        gene_df = pd.DataFrame()
        eqtl_df = pd.DataFrame()
    else:
        snp_df = pd.concat(all_snps)
        gene_df = pd.concat(all_genes)
        eqtl_df = pd.concat(all_eqtls)
    logger.verbose = True
    logger.write('  Time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60))
    return snp_df, gene_df, eqtl_df

def find_gene_snps(
        inter_df,
        spatial,
        enzyme,
        snps
        ):
    db.dispose()
    eqtl_project_db.dispose()

    if spatial == 'hic':
        table = 'variant_lookup_{}'
    else:
        table = 'variant_lookup_pchic_{}'
    
    if enzyme in ['MboI', 'DpnII']:  # MboI and DpnII have the same restriction sites
        table = table.format('mboi')
    else:
        table = table.format(enzyme.lower())
    
    if spatial == 'hic':
        chrom = inter_df['fragment_chr'].drop_duplicates().tolist()[0]
        sql = '''SELECT * FROM {}  WHERE chrom = '{}' AND frag_id >= {} AND frag_id <= {}'''
        df = pd.DataFrame()
        con = eqtl_project_db.connect()
        res = con.execute(
            sql.format(
                table, chrom, inter_df['fragment'].min(), inter_df['fragment'].max())
        ).fetchall()
        con.close()
        if res:
            df = pd.DataFrame(res, columns=['frag_id', 'chrom', 'id'])
        inter_df = inter_df.rename(columns={'chrom': 'gene_chr'})
        df = inter_df.merge(
            df, how='inner', left_on=['fragment'], right_on=['frag_id'])
        df['id'] = df['id'].astype('Int64')
        df = df[df['frag_id'].notnull()]
        df = df.drop_duplicates()
        snp_df = find_snp_by_id(spatial, df, eqtl_project_db)
        if snp_df.empty:
            return
        snp_df = snp_df.merge(df, how='inner', on=['id', 'chrom'])
        snp_df['enzyme'] = enzyme
        snp_df = snp_df.rename(columns={'rsid': 'snp'})
        snps.append(snp_df.drop_duplicates())
    else:
        snp_frag = inter_df['inter_frag'].drop_duplicates().tolist()
        if len(snp_frag) > 1:
            snp_frag = tuple(snp_frag)
            sql = '''SELECT * FROM {} WHERE frag_id IN {}'''
        else:
            snp_frag = (snp_frag[0])
            sql = '''SELECT * FROM {} WHERE frag_id = {}'''
        df = pd.DataFrame()
        con = eqtl_project_db.connect()
        res = con.execute(
                sql.format(
                    table, snp_frag)
                ).fetchall()
        con.close()
        if res:
            df = pd.DataFrame(res, columns=['frag_id', 'chrom', 'id'])
        inter_df = inter_df.rename(columns={'chrom': 'gene_chr'})
        df = inter_df.merge(
                df, how='inner', left_on=['inter_frag'], right_on=['frag_id'])
        df['id'] = df['id'].astype('Int64')
        df = df[df['frag_id'].notnull()]
        df = df.drop_duplicates()
        snp_df = find_snp_by_id(spatial, df, eqtl_project_db)
        if snp_df.empty:
            return
        snp_df = snp_df.merge(df, how='inner', on=['id', 'chrom'])
        snp_df['enzyme'] = enzyme
        snp_df = snp_df.rename(columns={'rsid': 'snp'})
        snps.append(snp_df.drop_duplicates())
         
def find_snp_by_id(spatial, df, db):
    df = df.sort_values(by=['id'])
    if spatial == 'hic':
        chunksize = eqtls.calc_chunksize(
            df['id'].tolist(), 2000000)
    else:
        chunksize = 100000
    chunks = [df[i:i+chunksize] for i in range(0, len(df), chunksize)]
    snp_df = []
    sql = '''SELECT rsid, variant_id, chrom, locus, id FROM variants WHERE id >= {}
    and id <= {}'''
    con = db.connect()
    for chunk in chunks:
        res = con.execute(sql.format(chunk['id'].min(), chunk['id'].max())).fetchall()
        if not res:
            continue
        res = pd.DataFrame(res, columns=[
            'rsid', 'variant_id', 'chrom', 'locus', 'id']).drop_duplicates()
        snp_df.append(res[res['id'].isin(chunk['id'])])
    con.close()
    if len(snp_df) == 0:
        return pd.DataFrame()
    return pd.concat(snp_df)

def find_snp_by_variant_id(df, db):
    snp_df = []
    sql = '''SELECT rsid, variant_id, chrom, locus, id FROM variants WHERE 
    variant_id >= '{}' and variant_id <= '{}' '''
    with db.connect() as con:
        snp_df = pd.read_sql(
            sql.format(df['variant_id'].min(), df['variant_id'].max()), con=con)
    if len(snp_df) == 0:
        return pd.DataFrame()
    snp_df = snp_df[snp_df['variant_id'].isin(df['variant_id'])]
    omitted_snps =  df[
        ~df['variant_id'].isin(snp_df['variant_id'])
    ]['variant_id'].drop_duplicates().tolist()
    return snp_df, omitted_snps

def filter_snp_fragments_hic(
        snp_df,
        logger):
    ''' Filter snp-fragment interactions '''
    logger.write('  * Filtering gene-SNP interactions...')
    snp_df['interactions'] = snp_df.groupby(
        ['variant_id', 'gene', 'cell_line', 'enzyme'])[
            'gene_fragment'].transform('count')
    snp_df['replicates'] = snp_df.groupby(
        ['variant_id', 'gene', 'cell_line', 'enzyme'])[
            'replicate'].transform('count')
    snp_df = snp_df.drop(columns=['replicate', 'gene_fragment'])
    snp_df = snp_df.drop_duplicates()
    snp_df['sum_interactions'] = snp_df.groupby(
        ['variant_id', 'gene'])[
            'interactions'].transform('sum')
    snp_df['sum_replicates'] = snp_df.groupby(
        ['variant_id', 'gene'])[
            'replicates'].transform('sum')
    snp_df['sum_cell_lines'] = snp_df.groupby(
        ['variant_id', 'gene'])['cell_line'].transform('count')
    condition = (
        (snp_df['interactions'] / snp_df['cell_line_replicates'] <= 1) &
        (snp_df['sum_replicates'] < 2) &
        (snp_df['sum_cell_lines'] < 2))
    gene_df = snp_df[~condition]
    snp_df = gene_df[
        ['id', 'snp', 'chrom', 'locus', 'variant_id', 'fragment', 'enzyme']
    ]
    snp_df.drop_duplicates(inplace=True)
    cols = ['snp', 'chrom', 'locus', 'variant_id',
            'gene', 'gencode_id', 'gene_chr', 'gene_start', 'gene_end',
            'interactions', 'replicates', 'enzyme',
            'cell_line', 'cell_line_replicates', 'sum_interactions',
            'sum_replicates', 'sum_cell_lines']
    return gene_df[cols], snp_df

def filter_snp_fragments_pchic(
        snp_df,
        logger):
    ''' Calculating N_reads and chicago_scores '''
    logger.write('  * Calculating cumulative scores for gene promoter-SNP interactions...')
    snp_df = snp_df.drop_duplicates()
    snp_df = snp_df.drop(columns=['replicate','id']).drop_duplicates()
    snp_df['N_reads'] = snp_df.groupby(
            ['variant_id', 'gencode_id', 'inter_frag',
            'gene_fragment', 'cell_line'])[
            'n_reads'].transform('sum')
    snp_df['Score'] = snp_df.groupby(['variant_id', 'gencode_id', 'inter_frag',
            'gene_fragment', 'cell_line','N_reads'])[
            'score'].transform('mean').round(2)
    gene_df = snp_df[['snp', 'chrom', 'locus', 'variant_id',
        'gene', 'gencode_id', 'gene_chr', 'gene_start', 'gene_end',
        'enzyme', 'cell_line', 'N_reads', 'Score']].drop_duplicates()
    snp_df = snp_df[
            ['snp', 'chrom', 'locus', 'variant_id','frag_id','enzyme']]
    snp_df.drop_duplicates(inplace=True)
    cols = ['snp', 'chrom', 'locus', 'variant_id',
            'gene', 'gencode_id', 'gene_chr', 'gene_start', 'gene_end',
            'enzyme', 'cell_line', 'N_reads', 'Score']
    return gene_df[cols], snp_df

def process_rs_df_whole(rs_df, db):
    rsid_tuple = ()
    if len(rs_df) == 1:
        rsid_tuple = str(tuple(rs_df['snp'])).replace(',', '')
    else:
        rsid_tuple = str(tuple(rs_df['snp']))
    sql = 'SELECT id, rsid, chrom, locus, variant_id FROM variants WHERE rsid IN {}'
    with db.connect() as con:
        return pd.read_sql_query(sql.format(rsid_tuple), con=con)


def process_rs_df(rs_df, db):
    rsid_tuple = ()
    df = []
    sql = '''SELECT id, rsid, chrom, locus, variant_id FROM variants WHERE rsid = '{}' '''
    with db.connect() as con:
        for idx, row in rs_df.iterrows():
            df.append(pd.read_sql_query(sql.format(row['snp']), con=con))
    df = pd.concat(df)
    return df


def process_position(position, db):
    chrom = ''
    locus = None
    temp_df = None
    try:
        chrom = position.split(':')[0].strip()
        locus = position.split(':')[1].strip()
        if '-' in locus:
            start = locus.split('-')[0].strip()
            end = locus.split('-')[1].strip()
            sql = '''SELECT id, rsid, chrom, locus, variant_id FROM variants
            WHERE chrom ='{}' AND locus >= {} AND locus <= {}'''
            with db.connect() as con:
                temp_df = pd.read_sql_query(sql.format(chrom, start, end), con=con)
        else:
            sql = '''SELECT id, rsid, chrom, locus, variant_id FROM variants
            WHERE chrom ='{}' AND locus = {}'''
            with db.connect() as con:
                temp_df = pd.read_sql_query(sql.format(chrom, locus), con=con)
        if not temp_df.empty:
            return temp_df
        else:
            return pd.DataFrame.from_dict(
                {'id': [np.nan], 'rsid': [position], 'chr': [np.nan],
                 'locus': [np.nan], 'variant_id': [np.nan]})
    except IndexError:
        return pd.DataFrame.from_dict({'id': [np.nan], 'rsid': [position], 'chr': [np.nan],
                                       'locus': [np.nan], 'variant_id': [np.nan]})


def get_snp_fragments_whole(snp_df, restriction_enzymes, db):
    snp_df = snp_df.sort_values(by=['id'])
    fragment_df = []
    chunksize = 1000
    chunks = [snp_df[i:i+chunksize]
              for i in range(0, snp_df.shape[0], chunksize)]
    for enzyme in restriction_enzymes:
        table = 'variant_lookup_{}'
        if enzyme in ['MboI', 'DpnII']:  # MboI and DpnII have the same restriction sites
            table = table.format('mboi')
        else:
            table = table.format(enzyme.lower())
        with db.connect() as con:
            for chunk in chunks:
                sql = '''SELECT * FROM {} WHERE id >= {} AND id <= {}'''
                df = pd.read_sql(sql.format(
                    table, chunk['id'].min(), chunk['id'].max()), con=con)
                df['enzyme'] = enzyme
                fragment_df.append(df.drop(columns=['chrom']))
    fragment_df = pd.concat(fragment_df)
    snp_df = pd.merge(snp_df, fragment_df, how='inner', on=['id'])
    db.dispose()
    return snp_df


def get_snp_fragments(spatial, snp_df, restriction_enzymes, db):
    snp_df = snp_df.sort_values(by=['id'])
    fragment_df = []
    chunksize = 1000
    chunks = [snp_df[i:i+chunksize]
              for i in range(0, snp_df.shape[0], chunksize)]
    for enzyme in restriction_enzymes:
        if spatial == 'pchic':
            table = 'variant_lookup_pchic_{}'
        else:
            table = 'variant_lookup_{}'
        df = []
        if enzyme in ['MboI', 'DpnII']:  # MboI and DpnII have the same restriction sites
            table = table.format('mboi')
        else:
            table = table.format(enzyme.lower())
        with db.connect() as con:
            for idx, row in snp_df.iterrows():
                sql = '''SELECT * FROM {} WHERE id = {}; '''
                df.append(pd.read_sql(sql.format(table, row['id']), con=con))
        df = pd.concat(df)
        df['enzyme'] = enzyme
        fragment_df.append(df.drop(columns=['chrom']))
    fragment_df = pd.concat(fragment_df)
    snp_df = pd.merge(snp_df, fragment_df, how='inner', on=['id'])
    db.dispose()
    return snp_df


def get_stdin_snp_position(inp, db):
    snp_df = []
    omitted_snps = []
    if inp.startswith('rs'):
        sql = '''SELECT id, rsid, chrom, locus, variant_id FROM variants WHERE rsid ='{}' '''
        with db.connect() as con:
            temp_df = pd.read_sql_query(sql.format(inp), con=con)
            if not temp_df.empty:
                snp_df.append(temp_df)
            else:
                omitted_snps.append(inp)
    elif inp.startswith('chr'):
        temp_df = process_position(inp, db)
        snp_df.append(temp_df[~temp_df['id'].isnull()])
        not_found = temp_df[temp_df['id'].isnull()]
        omitted_snps += not_found.rsid.tolist()
    return snp_df, omitted_snps


def get_file_snp_position(df, db):
    snp_df = []
    omitted_snps = []
    df = df[~df['snp'].isna()]
    rs_df = df[df.snp.str.startswith('rs')]
    chr_df = df[df.snp.str.startswith('chr')]
    if not rs_df.empty:
        temp_df = process_rs_df(rs_df, db)
        snp_df.append(temp_df)
        not_found = rs_df[~rs_df['snp'].isin(temp_df['rsid'])]
        omitted_snps += not_found.snp.tolist()
    if not chr_df.empty:
        temp_df = []
        for i, row in chr_df.iterrows():
            temp_df.append(process_position(row['snp'], db))
        temp_df = pd.concat(temp_df)
        snp_df.append(temp_df[~temp_df['id'].isnull()])
        not_found = temp_df[temp_df['id'].isnull()]
        omitted_snps += not_found.rsid.tolist()
    return snp_df, omitted_snps


def check_rs_merged(snps, rs_merge_arch_fp):
    '''
    Check human_9606_b151_GRCh38p7_RsMergeArch.bcp to see if missing rsID
    is merged into another rsID
    '''
    df = pd.DataFrame(
        [snp for snp in snps],  # if snp.startswith('rs')
        columns=['snp'])
    df['no_rs'] = df['snp'].str[2:]
    archive = pd.read_csv(
        rs_merge_arch_fp, sep='\t', compression='gzip',
        names=['new_id', 'old_id'],
        dtype={'new_id': 'object', 'old_id': 'object'})
    archive = archive[archive['new_id'].isin(df['no_rs']) |
                      archive['old_id'].isin(df['no_rs'])]
    new_df = []
    new_archive = df.merge(
        archive, how='inner', left_on='no_rs', right_on='new_id')
    if not new_archive.empty:
        new_archive = new_archive[~new_archive['new_id'].isnull()].drop(
            columns=['no_rs', 'new_id'])
        new_archive.columns = ['old_rs', 'snp']
        new_df.append(new_archive)
    old_archive = df.merge(
        archive, how='inner', left_on='no_rs', right_on='old_id')
    if not old_archive.empty:
        old_archive = old_archive[~old_archive['old_id'].isnull()].drop(
            columns=['no_rs', 'old_id'])
        old_archive.columns = ['old_rs', 'snp']
        new_df.append(old_archive)
    if len(new_df) == 0:
        return pd.DataFrame(), []
    new_df = pd.concat(new_df)
    new_df['snp'] = 'rs' + new_df['snp'].astype(str)
    omitted_snps = df[~df['snp'].isin(new_df['old_rs'])]['snp'].tolist()
    return new_df.reset_index(drop=True), omitted_snps


def process_snp_input(inputs, output_dir, db,
                      rs_merge_arch_fp, logger):
    snp_df = []
    omitted_snps = []
    len_input_snps = 0
    stdin_snps = []
    merged_snps = pd.DataFrame()
    for inp in inputs:
        if os.path.isfile(inp):
            df = []
            with open(inp) as f:
                for line in f:
                    if not line.strip() == '':
                        df.append([line.strip().lower()])
            df = pd.DataFrame(df, columns=['snp']).drop_duplicates() #Handle duplicated SNPs in the input file
            len_input_snps += len(df)
            temp_snp_df, temp_omitted_snps = get_file_snp_position(df, db)
            snp_df += temp_snp_df
            omitted_snps += temp_omitted_snps
        else:
            inp = inp.split(',')[0]
            if inp in stdin_snps:
                continue
            stdin_snps.append(inp)
            len_input_snps += 1
            temp_snp_df, temp_omitted_snps = get_stdin_snp_position(inp, db)
            snp_df += temp_snp_df
            omitted_snps += temp_omitted_snps
    _omitted_snps = []
    if len(omitted_snps) > 0:
        merged_snps, _omitted_snps = check_rs_merged(
            omitted_snps, rs_merge_arch_fp)
        if not merged_snps.empty:
            msg = '  * Warning: {} SNPs have their rsIDs merged'
            temp_df = process_rs_df(merged_snps[['snp']], db)
            if not temp_df.empty:
                in_list = len(temp_df[temp_df['variant_id'].isin(
                    pd.concat(snp_df)['variant_id'])])
                if in_list > 0:
                    msg = msg + \
                        ', {} of which are already on the input list'.format(
                            in_list)
                    len_input_snps -= in_list
                snp_df.append(temp_df)
                temp_merged_snps = merged_snps[merged_snps['snp'].isin(
                    temp_df['rsid'])]
                temp_omitted_snps = merged_snps[~merged_snps['old_rs'].isin(
                    temp_merged_snps['old_rs'])]['old_rs'].drop_duplicates()
                _omitted_snps += temp_omitted_snps.tolist()
                logger.write(msg.format(
                    len(temp_merged_snps)) + '. See snps_merged.txt for details.')
                merged_snps = temp_merged_snps.rename(
                    columns={'old_rs': 'input_rsid', 'snp': 'merged_rsid'})
            else:
                _omitted_snps = omitted_snps
        else:
            _omitted_snps = omitted_snps
        if len(_omitted_snps) > 0:
            msg = '  * Warning: {} SNPs not found in eQTL database.'.format(
                len(_omitted_snps))
            logger.write(msg + ' See snps_removed.txt for details.')
    if len(snp_df) == 0:
        snp_df = pd.DataFrame()
    else:
        snp_df = pd.concat(snp_df).drop_duplicates()
    return snp_df, len_input_snps, _omitted_snps, merged_snps

def get_snps_within_gene(gene, db):
    sql = '''SELECT id, rsid, chrom, locus, variant_id FROM variants 
    WHERE chrom = '{}' AND locus >='{}' AND locus <='{}' '''
    with db.connect() as con:
        df = pd.read_sql_query(
            sql.format(
                gene.loc[0, 'chrom'],
                gene.loc[0, 'start'],
                gene.loc[0, 'end']),
            con=con)
        return df['rsid'].drop_duplicates().tolist()

    
def get_snp(spatial,
            inputs,
            gene_out,
            hic_df,
            pchic_df,
            output_dir,
            db,
            rs_merge_arch_fp,
            logger,
            suppress_intermediate_files=False
            ):
    """Retrieve SNP position and restriction fragments.
    Args:
    inputs: File(s) (or stdin) containing SNP rsIDs or genomic positions in bed format (chr:start-end)>
    restriction_enzymes: a list of restriction enzymes with which query Hi-C/PCHi-C libraries are prepared
    output_dir: User-specified directory for results. Defaults to inputs directory.
    postgres_url: path to codes3d_common database
    suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
    A pandas dataframe
    with the df columns:
    1. snp #rsID
    2. SNP chromosome
    4. Fragment ID
    5. Fragment restriction enzyme
    If suppress_intermediate_files=False, write dataframe to snps.txt
    """
    logger.write('Processing SNP input...')
    start_time = time.time()
    if spatial == 'hic':
        restriction_enzymes = hic_df['enzyme'].drop_duplicates().tolist()
    else:
        restriction_enzymes = pchic_df['enzyme'].drop_duplicates().tolist()
    if gene_out is not None:
        if not gene_out.empty:
            inputs = get_snps_within_gene(gene_out, db)
    snp_df, len_input_snps, omitted_snps, merged_snps = process_snp_input(
        inputs, output_dir, db, rs_merge_arch_fp,  logger)
    if snp_df.empty:
        logger.write('We could not find your SNPs in our databases.')
        sys.exit()
    snp_df[['id', 'locus']] = snp_df[['id', 'locus']].astype(int)
    snp_df = get_snp_fragments(spatial, snp_df, restriction_enzymes, db)
    snp_df = snp_df.rename(columns={'rsid': 'snp', 'frag_id': 'fragment'})
    if not suppress_intermediate_files:
        if not merged_snps.empty:
            merged_snps.to_csv(os.path.join(output_dir, 'snps_merged.txt'),
                               sep='\t', index=False)
        if len(omitted_snps) > 0:
            pd.DataFrame(omitted_snps).to_csv(
                os.path.join(output_dir, 'snps_removed.txt'),
                sep='\t', index=False, header=False)
    logger.write('  * {} out of {} SNPs passed parsing criteria.'.format(
        len(snp_df['snp'].drop_duplicates()), len_input_snps))
    logger.write('  Time elasped: {:.2f} mins'.format(
        (time.time()-start_time)/60))
    db.dispose()
    return snp_df

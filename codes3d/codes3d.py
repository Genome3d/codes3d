#!/usr/bin/env python

import os
import sys
import pandas as pd
import argparse
import configparser
import multiprocessing
import time
import datetime
from sqlalchemy import create_engine
from sqlalchemy.pool import NullPool
import tqdm
import statsmodels.stats.multitest as multitest
import snps
import genes
import interactions
import summary
import eqtls
import aFC

def parse_tissues(user_tissues, match_tissues, eqtl_project, db):
    if eqtl_project:
        sql = '''SELECT * FROM meta_eqtls WHERE project = '{}' '''.format(
            eqtl_project)
    else:
        sql = '''SELECT * FROM meta_eqtls'''
    df = pd.DataFrame()
    with db.connect() as con:
        df = pd.read_sql(sql, con=con)
    db.dispose()
    tissues = []
    if match_tissues:
        user_tissues = match_tissues[0]
    if user_tissues:
        matched_df = []
        matched_tissues = []
        to_omit = []
        not_tissues = []
        for u_tissue in user_tissues:
            u_df = df[
                (df['name'] == u_tissue) |
                (df['tags'].str.contains(
                    r'\b{}\b'.format(u_tissue), case=False))
            ]
            if u_df.empty:
                if u_tissue.startswith('-'):
                    to_omit.append(u_tissue)
                else:
                    not_tissues.append(u_tissue)
            else:
                matched_df.append(u_df)
                matched_tissues.append(u_tissue)
        error_msg = 'Program aborting:\n\t{}\nnot found in database.'
        if (len(matched_df) == 0 or len(not_tissues) > 0) and len(to_omit) == 0:
            print(error_msg.format('\n\t'.join(not_tissues)))
            print('\nPlease use one of the following. ' +
                  'Tissue names are case sensitive:')
            list_eqtl_tissues(db)
            sys.exit()
        user_df = pd.DataFrame()
        if len(to_omit) > 0 and len(matched_tissues) == 0:
            user_df = df
        else:
            user_df = pd.concat(matched_df)
        if match_tissues:
            for i in range(len(matched_tissues)):
                user_df = user_df[
                    user_df['tags'].str.contains(
                        r'\b{}\b'.format(matched_tissues[i]), case=False)]
            user_df = user_df.drop_duplicates()
            for i in range(len(to_omit)):
                user_df = user_df[
                    ~user_df['tags'].str.contains(
                        r'\b{}\b'.format(to_omit[i][1:]), case=False)]
        if len(user_df['project'].drop_duplicates()) > 1 and not eqtl_project:
            # Ensure tissues are from same eQTL project
            print('FATAL: eQTL tissues are from different projects. ',
                  'Add another tag to fine-tune match',
                  'or use \'--eqtl-project\' to specify project.')
            print(user_df[['name', 'project']].to_string(index=False))
            sys.exit()
        tissues = user_df[['name', 'project']]
    else:  # Use GTEx database as default
        tissues = df[df['project'] == 'GTEx'][[
            'name', 'project']]
    return tissues

def parse_hic(
        match_tissues,
        include_cell_line,
        exclude_cell_line,
        restriction_enzymes,
        db,
        pchic=False):
    ''' user parameters -r, -n and -x.

    Args:
        restriction_enzymes: space-delimited list of restriction enzymes from
        user. Limits program to Hic libraries restricted by specified enzyme.
        include_cell_line: space-delimited list of cell_lines from -n.
        exclude_cell_line: space-delimited list of cell_lines from -x
    Returns:
        hic_df: a dataframe columns(library, enzyme, rep_count)
        pchic_df: a dataframe columns(library, enzyme, rep_count)
    '''
    if not args.pchic:
        sql = '''SELECT library, tags, enzyme, rep_count FROM meta_hic'''
    else:
        sql = '''SELECT library, tags, enzyme, rep_count FROM meta_pchic'''
    
    df = pd.DataFrame()
    with db.connect() as con:
        df = pd.read_sql(sql, con=con)
    db.dispose()
    hic_df = pd.DataFrame()
    if match_tissues:
        matched_df = []
        matched_tissues = []
        to_omit = []
        not_tissues = []
        for u_tissue in match_tissues[0]:
            u_df = df[
                    (df['library'] == u_tissue) |
                    (df['tags'].str.contains(
                        r'\b{}\b'.format(u_tissue), case=False))
            ]
            #u_df = df[
            #        df['library'].str.contains(u_tissue,case=False) | 
            #        df['tags'].str.contains(u_tissue,case=False)]
            if u_df.empty:
                if u_tissue.startswith(u_tissue):
                    to_omit.append(u_tissue)
                else:
                    not_tissues.append(u_tissue)
            else:
                matched_df.append(u_df)
                matched_tissues.append(u_tissue)
        error_msg = 'Program aborting:\n\t{}\ndid not match any Hi-C library.'
        if (len(matched_df) == 0 or len(not_tissues) > 0) and len(to_omit) == 0:
            print(error_msg.format('\n\t'.join(not_tissues)))
            print(('Use -t and -n to include specific eQTL tissues'
                ' and Hi-C/PCHi-C libraries. Library names are case sensitive:'))
            sys.exit('\n\t{}'.format('\n\t'.join(df['library'].tolist())))
        if len(matched_df) == 0 and len(to_omit) > 0:
            hic_df = df
        else:
            hic_df = pd.concat(matched_df)
        if match_tissues:
            for i in range(len(matched_tissues)):
                hic_df = hic_df[
                        hic_df['tags'].str.contains(
                            r'\b{}\b'.format(matched_tissues[i]), case=False)]
            for i in range(len(to_omit)):
                hic_df = hic_df[
                        ~hic_df['tags'].str.contains(
                        r'\b{}\b'.format(to_omit[i][1:]), case=False)]
            hic_df = hic_df.drop_duplicates()
    elif include_cell_line and len(include_cell_line) > 0:
        validate_input(include_cell_line, df['library'])
        hic_df = df[df['library'].str.upper().isin(
            [c.upper() for c in include_cell_line])]
    elif exclude_cell_line and len(exclude_cell_line) > 0:
        validate_input(exclude_cell_line, df['library'])
        hic_df = df[~df['library'].str.upper().isin(
            [c.upper() for c in exclude_cell_line])]
    else:
        hic_df = df
    if restriction_enzymes and len(restriction_enzymes) > 0:
        validate_input(restriction_enzymes, df['enzyme'].drop_duplicates())
        hic_df = hic_df[hic_df['enzyme'].str.upper().isin(
            [c.upper() for c in restriction_enzymes])]
    if not hic_df.empty:
        return hic_df.drop_duplicates()
    else:
        return df.drop_duplicates()

def validate_input(input_params, params):
    input_params_upper = [c.upper() for c in input_params]
    not_found = set(input_params_upper).difference(
        set(params.str.upper().tolist()))
    if len(not_found) > 0:
        print('FATAL: The following parameters are not recognized:')
        for param in not_found:
            print('\t{}'.format(input_params[input_params_upper.index(param)]))
        sys.exit(
            'Please ensure parameter is from the following:\n\t{}'.format(
                '\n\t'.join(params.tolist())))


def calc_afc(eqtl_df, genotypes_fp, expression_dir, covariates_dir,
             eqtl_project, output_dir, fdr_threshold,  bootstrap, num_processes):
    
    if 'adj_pval' in eqtl_df.columns: # Exclude non-significant eQTLs.
        eqtl_df = eqtl_df[eqtl_df['adj_pval'] <= fdr_threshold]
    if eqtl_df.empty:
        print('Warning: No significant eQTL associations found.\nExiting.')
        sys.exit()
    eqtl_df = aFC.main(
        eqtl_df, genotypes_fp, expression_dir, covariates_dir,
        eqtl_project, output_dir, bootstrap, num_processes)

    return eqtl_df


def list_eqtl_databases(db):
    sql = '''SELECT * FROM meta_eqtls'''
    with db.connect() as con:
        df = pd.read_sql(sql, con=con)
        df = df[['project']].drop_duplicates().reset_index()
        for idx, row in df.iterrows():
            print('{}\t{}'.format(idx + 1, row['project']))
    db.dispose()


def list_eqtl_tissues(db):
    sql = '''SELECT * FROM meta_eqtls'''
    with db.connect() as con:
        df = pd.read_sql(sql, con=db)
        for idx, row in df.iterrows():
            print('{}\t{}'.format(idx + 1, row['name']))
    db.dispose()

def list_hic_libraries(db):
    sql = '''SELECT library, tissue FROM meta_hic'''
    with db.connect() as con:
        df = pd.read_sql(sql, con=db)
        for idx, row in df.drop_duplicates().iterrows():
            print('{}\t{}'.format(idx + 1, row['library']))
    db.dispose()

def list_pchic_libraries(db):
    sql = '''SELECT library, tissue FROM meta_hic'''
    with db.connect() as con:
        df = pd.read_sql(sql, con=db)
        for idx, row in df.drop_duplicates().iterrows():
            print('{}\t{}'.format(idx + 1, row['library']))
    db.dispose()

def list_enzymes(db):
    sql = '''SELECT DISTINCT enzyme FROM meta_hic'''
    with db.connect() as con:
        df = pd.read_sql(sql, con=db)
        for idx, row in df.drop_duplicates().iterrows():
            print('{}\t{}'.format(idx + 1, row['enzyme']))
    db.dispose()

def list_enzymes_pchic(db):
    sql = '''SELECT DISTINCT enzyme FROM meta_pchic'''
    with db.connect() as con:
        df = pd.read_sql(sql, con=db) 
        for idx, row in df.drop_duplicates().iterrows():
            print('{}\t{}'.format(idx + 1, row['enzyme']))
    db.dispose()

def list_tissue_tags(db):
    sql = '''SELECT tags FROM {}'''
    hic_df = pd.DataFrame()
    eqtl_df = pd.DataFrame()
    with db.connect() as con:
        hic_df = pd.read_sql(sql.format('meta_hic'), con=con)
        eqtl_df = pd.read_sql(sql.format('meta_eqtls'), con=con)
    db.dispose()
    hic_tags = []
    eqtl_tags = []
    for idx, row in hic_df.iterrows():
        hic_tags += row['tags'].split(', ')
    for idx, row in eqtl_df.iterrows():
        eqtl_tags += row['tags'].split(', ')
    tags = sorted(list(set(hic_tags).intersection(set(eqtl_tags))))
    print(pd.DataFrame(tags, columns=['tags']
                      ).to_string(index=False, header=None))

def list_tissue_tags_pchic(db):
    sql = '''SELECT tags FROM {}'''
    pchic_df = pd.DataFrame()
    eqtl_df = pd.DataFrame()
    with db.connect() as con:
        pchic_df = pd.read_sql(sql.format('meta_hic'), con=con)
        eqtl_df = pd.read_sql(sql.format('meta_eqtls'), con=con)
    db.dispose()
    pchic_tags = []
    eqtl_tags = []
    for idx, row in pchic_df.iterrows():
        pchic_tags += row['tags'].split(', ')
    for idx, row in eqtl_df.iterrows():
        eqtl_tags += row['tags'].split(', ')
    tags = sorted(list(set(pchic_tags).intersection(set(eqtl_tags))))
    print(pd.DataFrame(tags, columns=['tags']
                       ).to_string(index=False, header=None))


def parse_intermediate_files(inter_files, output_dir, file_type, multi_test=['all']):
    df = []
    for inter_file in tqdm.tqdm(
            inter_files,
            desc='Parsing {} files'.format(file_type),
            unit=' file'):
        for chunk_df in pd.read_csv(inter_file, sep='\t', chunksize=1000000):
            df.append(chunk_df)
    df = pd.concat(df).drop_duplicates()
    if file_type == 'eqtls':
        df = multi_test_correction(df, multi_test)
    return df

def multi_test_correction(eqtl_df, multi_test):
    print('  * Adjusting eQTL pvalues...')
    cols = ['sid', 'pid', 'sid_chr', 'sid_pos',
            'adj_pval', 'pval', 'b', 'b_se', 'maf', 'tissue']
    if 'snp' == multi_test.lower():
        adj_pval_df = []
        for _, df in eqtl_df.groupby(['sid', 'tissue']):
            df['adj_pval'] = correct_pvals(df['pval'])
            adj_pval_df.append(df)
        adj_pval_df = pd.concat(adj_pval_df)
        return adj_pval_df[cols]
    if 'tissue' == multi_test.lower():
        adj_pval_df = []
        for _, df in eqtl_df.groupby(['tissue']):
            df['adj_pval'] = correct_pvals(df['pval'])
            adj_pval_df.append(df)
        adj_pval_df = pd.concat(adj_pval_df)
        return adj_pval_df[cols]
    else:
        eqtl_df['adj_pval'] = correct_pvals(eqtl_df['pval'])
        return eqtl_df[cols]

def correct_pvals(pval):
    return multitest.multipletests(pval, method='fdr_bh')[1]


def map_gtex_cis_eqtls(
        snp_df,
        tissues,
        C,
        args,
        eqtl_project_db,
        logger):
    eqtl_df = eqtls.map_eqtls_non_spatial(
        snp_df,
        tissues,
        C.eqtl_data_dir,
        args.num_processes,
        eqtl_project_db,
        logger)
    gene_df = genes.get_gene_by_gencode(
        eqtl_df.rename(columns={'gene_id': 'gene'}),
        commons_db)[0]
    gene_df = pd.concat(gene_df).rename(
        columns = {
            'name': 'gene',
            'chrom': 'gene_chrom',
            'start': 'gene_start',
            'end': 'gene_end'}).drop(
                columns = ['id'])
    eqtl_df = eqtl_df.merge(
        gene_df, how = 'inner',
        left_on = ['gene_id'], right_on = ['gencode_id']
    ).drop(columns = ['gene_id']).drop_duplicates()
    cols = ['snp', 'variant_id', 'gene', 'gencode_id',
            'pval_nominal', 'slope', 'slope_se', 'pval_nominal_threshold','min_pval_nominal', 'pval_beta',
            'tss_distance', 'maf', 'gene_chrom', 'gene_start', 'gene_end']

    if not args.no_afc:
        afc_start_time = time.time()
        eqtl_df = calc_afc(
            eqtl_df,
            genotypes_fp,
            expression_dir,
            covariates_dir,
            eqtl_project,
            args.output_dir,
            args.fdr_threshold,
            args.afc_bootstrap,
            args.num_processes)

    print(eqtl_df)
    eqtl_df[cols].to_csv(os.path.join(args.output_dir, 'non_spatial_eqtls.txt'), sep='\t', index=False)
    msg = 'Done.\nTotal time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60)
    logger.write(msg)
    
    
    sys.exit()
    
def map_non_spatial_eqtls(
        snp_df,
        gene_df,
        tissues,
        hic_df,
        C,
        args,
        genotypes_fp,
        eqtl_project_db,
        commons_db,
        covariates_dir,
        expression_dir,
        logger,
        pchic = False):
    gene_info_df = None
    gene_list = None
    if not snp_df.empty:
        if args.gene_list is None:
            print('No gene list given. Will map across all genes.')
        else:
            gene_info_df = genes.get_gene_info(
                    args.gene_list,
                    hic_df,
                    args.output_dir,
                    commons_db,
                    logger,
                    pchic,
                    args.suppress_intermediate_files)
            gene_info_df = gene_info_df.rename(columns={
                'name': 'gene',
                'chrom': 'gene_chrom',
                'start': 'gene_start',
                'end': 'gene_end'})
            gene_list = gene_info_df['gencode_id'].drop_duplicates().tolist()
    if not gene_df.empty:
        gene_info_df = gene_df
        gene_list = gene_info_df['gencode_id'].drop_duplicates().tolist()
    if args.gtex_cis:
        eqtl_df = eqtls.map_cis_eqtls_gtex(
            snp_df,
            gene_list,
            tissues,
            C.eqtl_data_dir,
            args.num_processes,
            eqtl_project_db,
            logger)
    else:
        eqtl_df = eqtls.map_eqtls_non_spatial(
            snp_df,
            gene_list,
            tissues,
            args.output_dir,
            args.maf_threshold,
            C,
            args.num_processes,
            genotypes_fp,
            eqtl_project_db,
            covariates_dir,
            expression_dir,
            logger)
    if gene_info_df is None:
        gene_info_df = genes.get_gene_by_gencode(
            eqtl_df.rename(columns={'phenotype_id': 'gene'}),
            commons_db)[0]
        gene_info_df = pd.concat(gene_info_df).rename(
            columns = {
                'name': 'gene',
                'chrom': 'gene_chrom',
                'start': 'gene_start',
                'end': 'gene_end'}).drop(
                    columns = ['id'])
    if snp_df.empty:
        snp_df = snps.find_snp_by_variant_id(
            eqtl_df[['variant_id']].drop_duplicates(), eqtl_project_db)[0]
        #snp_df = pd.concat(snp_df)
        snp_df = snp_df.rename(columns={'rsid': 'snp'})
    gene_info_df = gene_info_df.rename(columns={
        'name': 'gene',
        'chrom': 'gene_chrom',
        'start': 'gene_start',
        'end': 'gene_end'})
    eqtl_df = eqtl_df.merge(
        gene_info_df, how = 'inner',
        left_on = ['phenotype_id'], right_on = ['gencode_id']).drop_duplicates()
    snp_df = snp_df.rename(columns={
        'chrom': 'snp_chrom',
        'locus': 'snp_locus'})
    eqtl_df = eqtl_df.merge(snp_df, how = 'left',
                            on=['variant_id'])
    cols = ['snp', 'variant_id', 'gene', 'gencode_id', 'tissue',
            'pval', 'b', 'b_se',
            'snp_chrom', 'snp_locus', 'maf',
            'gene_chrom', 'gene_start', 'gene_end']
    eqtl_df = eqtl_df[cols].drop_duplicates()
    eqtl_project = tissues['project'].iloc[0]
    if not args.no_afc:
        afc_start_time = time.time()
        eqtl_df['sid'] = eqtl_df['variant_id']
        eqtl_df['sid_chr'] = eqtl_df['snp_chrom']
        eqtl_df['sid_pos'] = eqtl_df['snp_locus']
        eqtl_df['pid'] = eqtl_df['gencode_id']
        eqtl_df = calc_afc(
            eqtl_df,
            genotypes_fp,
            expression_dir,
            covariates_dir,
            eqtl_project,
            args.output_dir,
            args.fdr_threshold,
            args.afc_bootstrap,
            args.num_processes)
        cols += ['log2_aFC', 'log2_aFC_lower', 'log2_aFC_upper']
    fp = os.path.join(args.output_dir, 'non_spatial_eqtls.txt')
    eqtl_df[cols].to_csv(fp, sep='\t', index=False)
    logger.write(f'Output written to {fp}')
    msg = 'Done.\nTotal time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60)
    logger.write(msg)
    sys.exit()

def map_spatial_eqtls(
        snp_list,
        C,
        hic_df,
        args,
        logger,
        commons_db,
        tissues,
        genotypes_fp,
        eqtl_project_db,
        covariates_dir,
        expression_dir,
        pchic=False):
    gene_df = []
    eqtl_df = []
    batchsize = 2000
    snp_batches = [
        snp_list[i:i + batchsize] for i in range(0, len(snp_list),
                                                 batchsize)]
    for batch_num, snp_batch in enumerate(snp_batches):
        #batch_output_dir = args.output_dir
        if len(snp_batches) > 1:
            logger.write('SNP batch {} of {}'.format(
                batch_num+1, len(snp_batches)))
            #batch_output_dir = os.path.join(
            #    args.output_dir, str(batch_num))
            #os.makedirs(batch_output_dir, exist_ok=True)
        batch_snp_df = snp_df[snp_df['snp'].isin(snp_batch)]
        batch_interactions_df = interactions.find_interactions(
            batch_snp_df,
            'snp',
            C.lib_dir,
            hic_df,
            args.num_processes,
            logger,
            pchic)
        if pchic:
            batch_gene_df = genes.get_gene_by_fid(
                batch_snp_df,
                batch_interactions_df,
                commons_db,
                logger,
                pchic)
        else:
            batch_gene_df = genes.get_gene_by_id(
                batch_snp_df,
                batch_interactions_df,
                commons_db,
                logger)
        gene_df.append(batch_gene_df)
        if batch_gene_df.empty:
            continue
        batch_eqtl_df = eqtls.map_eqtls(
            batch_gene_df,
            tissues,
            args.output_dir,
            C,
            genotypes_fp,
            args.num_processes,
            eqtl_project_db,
            covariates_dir,
            expression_dir,
            args.pval_threshold,
            args.maf_threshold,
            args.fdr_threshold,
            logger)
        eqtl_df.append(batch_eqtl_df)
    gene_df = pd.concat(gene_df)
    eqtl_df = pd.concat(eqtl_df)
    if gene_df.empty or snp_df.empty or eqtl_df.empty:
        logger.write('''No eQTLs found.Cleaning up genes.txt and eqtls.txt.
        \nProgram exiting.''')
        sys.exit()
    if not args.suppress_intermediate_files:
        gene_df.to_csv(os.path.join(
            args.output_dir, 'genes.txt'), sep='\t', index=False)

    return gene_df, eqtl_df

class Logger(object):
    def __init__(self, logfile=None, verbose=True):
        """
        A simple logger.
        """
        self.console = sys.stdout
        self.verbose = verbose
        if logfile is not None:
            self.log = open(logfile, 'w')
        else:
            self.log = None

    def write(self, message):
        if self.verbose:
            self.console.write(message+'\n')
        if self.log is not None:
            self.log.write(message+'\n')
            self.log.flush()

    def verbose(self, verbose=True):
        self.verbose = verbose


class CODES3D:
    def __init__(self, config_fp):
        config = configparser.ConfigParser()
        config.read(config_fp)
        self.lib_dir = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults", "LIB_DIR"))
        self.hic_dir = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults", "HIC_DIR"))
        self.pchic_dir = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults", "PCHIC_DIR"))
        self.gene_bed_fp = os.path.join(
            os.path.dirname(__file__), config.get("Defaults", "GENE_BED_FP"))
        self.snp_bed_dir = os.path.join(
            os.path.dirname(__file__), config.get("Defaults", "SNP_BED_DIR"))
        self.gene_database_fp = os.path.join(
            os.path.dirname(__file__), config.get("Defaults", "GENE_DATABASE_FP"))
        self.eqtl_data_dir = os.path.join(
            os.path.dirname(__file__), config.get("Defaults", "EQTL_DATA_DIR"))
        self.rs_merge_arch_fp = os.path.join(
            os.path.dirname(__file__), config.get("Defaults", "RS_MERGE_ARCH"))
        self.host = config.get("postgresql", "host")
        self.plink = os.path.join(
            os.path.dirname(__file__), "plink2")
        self.commons_db = config.get("postgresql", "database")
        self.user = config.get("postgresql", "user")
        self.password = config.get("postgresql", "password")
        self.port = config.get("postgresql", "pgbouncer_port")
        self.port = 5432
        self.commons_db_url = 'postgresql://{}:{}@{}:{}/{}'.format(
            self.user, self.password, self.host, self.port, self.commons_db)
        self.eqtl_db_url = 'postgresql://{}:{}@{}:{}/{}'.format(
            self.user, self.password, self.host, self.port, 'eqtls_{}')
        self.hic_restriction_enzymes = [
            e.strip() for e in os.listdir(self.hic_dir)]
        self.pchic_restriction_enzymes = [
            e.strip() for e in os.listdir(self.pchic_dir)]

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''CoDeS3D maps gene regulatory landscape using chromatin 
        interaction and eQTL data.''')
    parser.add_argument(
        '-s', '--snp-input', nargs='+',
        help='''The dbSNP IDs or loci of SNPs of interest in the format
         \'chr<x>:<locus>\'.\n
        Use this flag to identify eGenes associated with the SNPs of interest.''')
    parser.add_argument(
        '-g', '--gene-input', nargs='+',
        help='''The symbols, Ensembl IDs or loci of genes interest in the format
         \'chr<x>:<start>-<end>\'.\n
        Use this flag to identify eQTLs associated with the gene of interest.''')
    parser.add_argument(
        '--snps-within-gene', nargs='+',
        help='''A gene symbol, Ensembl ID or location in the format
        \'chr<x>:<start>-<end>\'.\n
        Use this flag to identify eGenes associated with the SNPs located within 
        the gene of interest.''')
    parser.add_argument(
        '-o', '--output-dir',
        help='The directory in which to output results.')
    parser.add_argument(
        '--pchic', action='store_true', default=False,
        help='''Use this flag to use pchic datasets instead of hic.''') 
    parser.add_argument(
        '--multi-test', default = 'multi',
        help='''Options for BH multiple-testing: ['snp', 'tissue', 'multi'].
        'snp': corrects for genes associated with a given SNP in a given tissue.
        'tissue': corrects for all associations in a given tissue.
        'multi': corrects for all associations across all tissues tested.''')
    parser.add_argument(
        '--pval-threshold', default=1,
        help='Maximum p value for mapping eQTLs. Default: 1.')
    parser.add_argument(
        '-f', '--fdr-threshold', type=float, default=0.05,
        help='''The FDR threshold to consider an eQTL
         statistically significant (default: 0.05).''')
    parser.add_argument(
        '--maf-threshold', default=0.1,
        help='Minimum MAF for variants to include. Default: 0.1.')
    parser.add_argument(
        '-p', '--num-processes', type=int,
        default=int(multiprocessing.cpu_count()/4),
        help='Number of CPUs to use (default: half the number of CPUs).')
    parser.add_argument(
        '--no-afc', action='store_true', default=False,
        help='''Do not calculate allelic fold change (aFC).
        If true, eQTL beta (normalised effect size) is calculated instead.
        (default: False).''')
    parser.add_argument(
        '--afc-bootstrap', type=int, default=1000,
        help='Number of bootstrap for  aFC calculation (default: 1000).')
    parser.add_argument(
        '-n', '--include-cell-lines', nargs='+',
        help='''Space-separated list of cell lines to include
               (others will be ignored). NOTE: Mutually exclusive
               with EXCLUDE_CELL_LINES. --match-tissues takes precedence.''')
    parser.add_argument(
        '-x', '--exclude-cell-lines', nargs='+',
        help='''Space-separated list of cell lines to exclude 
         (others will be included). NOTE: Mutually exclusive
         with INCLUDE_CELL_LINES. --match-tissues and -n take precedence.''')
    parser.add_argument(
        '--list-hic-libraries', action='store_true', default=False,
        help='List available Hi-C libraries.')
    parser.add_argument(
        '--list-pchic-libraries', action='store_true', default=False,
        help='List available PCHi-C libraries.')
    parser.add_argument(
        '--match-tissues', action='append', nargs=argparse.REMAINDER, default=None,
        help='''Try to match eQTL and Hi-C tissue types using space-separated 
        tags. When using this, make sure that it is the last tag. 
        Note that tags are combined with the AND logic. 
        Prepend "-" to a tag if you want it to be excluded.
        Use `--list-tissues-tags' for possible tags.''')
    parser.add_argument(
        '--list-tissue-tags', action='store_true',
        help='''List tags to be used with `--match-tissues'.''')
    parser.add_argument(
        '-t', '--tissues', nargs='+',
        help='''Space-separated list of eQTL tissues to query.
        Note that tissues are case-sensitive and must be from the same eQTL projects.
        Default is all tissues from the GTEx project.
        Use 'codes3d.py --list-eqtl-tissues' for a list of installed tissues.''')
    parser.add_argument(
        '--eqtl-project', type=str, default=None,
        help='''The eQTL project to query. Default: GTEx.
        'use \'codes3d.py --list-eqtl-db\' to list available databases.''')
    parser.add_argument(
        '--list-eqtl-db', action='store_true',
        help='List available eQTL projects to query.')
    parser.add_argument(
        '-r', '--restriction-enzymes', nargs='+',
        help='''Space-separated list of restriction enzymes used in Hi-C data.'
         Use \'list-enzymes\' to see available enzymes.''')
    parser.add_argument(
        '--list-enzymes', action='store_true',
        help='List restriction enzymes used to prepare installed Hi-C libraries.')
    parser.add_argument(
        '--list-enzymes-pchic', action='store_true',
        help='List restriction enzymes used to prepare installed PCHi-C libraries.')
    parser.add_argument(
        '--list-eqtl-tissues', action='store_true',
        help='List available eQTL tissues to query.')
    parser.add_argument(
        '-c', '--config',
        default=os.path.join(os.path.dirname(
            __file__), '../docs/codes3d.conf'),
        help='The configuration file to use (default: docs/codes3d.conf).')
    parser.add_argument(
        '--output-format', type=str, default='long',
        help='''Determines columns of output.\n 
        \'short\': snp|gencode_id|gene|tissue|adj_pval|
        [log2_aFC|log2_aFC_lower|log2_aFC_upper] or [beta|beta_se]|
        maf|interaction_type|hic_score.\n 
        \'medium\' adds: eqtl_pval|snp_chr|snp_locus|ref|alt|gene_chr|
        gene_start|gene_end|distance.\n 
         \'long\' adds:cell_lines|cell_line_hic_scores|
        expression|max_expressed_tissue|max_expression|
        min_expressed_tissue|min_expression.\n(default: long).''')
    parser.add_argument(
        '--do-not-produce-summary', action='store_true', default=False,
        help='Do not produce final summary file, stop process after ' +
        'mapping eQTLs. Helpful if running batches (default: False).')
    parser.add_argument(
        '--suppress-intermediate-files', action='store_true', default=False,
        help='''Do not produce intermediate 
         files. These can be used to run the pipeline from
         an intermediate stage in the event of interruption
         (default: False).''')
    parser.add_argument(
        '--non-spatial', action='store_true', default=False,
        help='Map non-spatial eQTLs.')
    parser.add_argument(
        '--gene-list', nargs='+', default=None,
        help='List of genes for non-spatial eQTL mapping.')
    parser.add_argument(
        '--gtex-cis', action='store_true', default=False,
        help='''Retrieve spatially unconstrained cis-eQTLs as calculated in GTEx.
        To be used in combination with 'non-spatial'. ''')
    
    return parser.parse_args()

def validate_args(args, commons_db):
    if args.list_eqtl_db:
        list_eqtl_databases(commons_db)
        sys.exit()
    if args.list_eqtl_tissues:
        list_eqtl_tissues(commons_db)
        sys.exit()
    if args.list_hic_libraries:
        list_hic_libraries(commons_db)
        sys.exit()
    if args.list_pchic_libraries:
        list_pchic_libraries(commons_db)
        sys.exit()
    if args.list_enzymes:
        list_enzymes(commons_db)
        sys.exit()
    if args.list_enzymes_pchic:
        list_enzymes_pchic(commons_db)
        sys.exit()
    if args.list_tissue_tags:
        list_tissue_tags(commons_db)
        sys.exit()
    if not (args.snp_input or args.gene_input or args.snps_within_gene) or not args.output_dir:
        print(
            '''Missing --snp-input, --gene-input, --gene-out, or --output-dir  
            parameter(s).''')
        sys.exit('\tUse \'codes3d.py -h\' for more details.')
    if (args.snp_input and args.gene_input) or\
       (args.snp_input and args.snps_within_gene)or \
       (args.gene_input and args.snps_within_gene):
        sys.exit('''FATAL: Use only one of --snp-input, --gene-input, or --gene-out''')

def log_settings(args, logger):
    now = datetime.datetime.now()
    logger.write(f'{now.strftime("%d/%m/%Y %H:%M:%S")}')
    mode = 'SNP'
    if args.gene_input:
        mode = 'Gene'
    elif args.snps_within_gene:
        mode = 'SNPs-within-gene'
    multi_test = 'Multi-tissue'
    if args.multi_test.lower() == 'snp':
        multi_test = 'SNP'
    elif args.multi_test.lower() == 'tissue':
        multi_test = 'Tissue'
    effect_size = 'allelic fold change (aFC)' if not args.no_afc else 'beta (normalised)'
    logger.write('SETTINGS')        
    logger.write(f'Run mode:\t{mode}')
    logger.write(f'FDR correction on:\t{multi_test}')
    logger.write(f'FDR threshold:\t{args.fdr_threshold}')
    logger.write(f'MAF threshold:\t{args.maf_threshold}')
    logger.write(f'Effect size:\t{effect_size}')

    if not args.pchic and not args.non_spatial:
        logger.write('Using Hi-C datasets\n')
    else:
        if args.pchic and not args.non_spatial:
            logger.write('Using PCHi-C datasets\n')

    if not args.pchic and not args.match_tissues and not args.include_cell_lines and \
        not args.exclude_cell_lines:
        if not args.non_spatial:
            logger.write('Hi-C libraries:\tAll libraries in database')
    else:
        if not args.pchic and not  args.non_spatial:
            logger.write('Hi-C libraries:\t{}'.format(
                ', '.join(hic_df['library'].tolist())))


    if args.pchic and not args.match_tissues and not args.include_cell_lines and \
       not args.exclude_cell_lines:
        if not args.non_spatial:
            logger.write('PCHi-C libraries:\tAll libraries in database')
    else:
        if args.pchic and not args.non_spatial:
            logger.write('PCHi-C libraries:\t{}'.format(
                ', '.join(hic_df['library'].tolist())))


    if not args.tissues and not args.match_tissues:
        logger.write('eQTL tissues:\tAll tissues in database\n')
    else:
        logger.write('\neQTL tissues:\t{}\n'.format(
            ', '.join(tissues['name'].tolist())))
    if args.snp_input:
        logger.write(f'--snp-input:\t{", ".join(args.snp_input)}')
    if args.gene_input:
        logger.write(f'--gene-input:\t{", ".join(args.gene_input)}')
    if args.snps_within_gene:
        logger.write(f'--snps-within-gene:\t{args.snps_within_gene}')
    logger.write('\n')


        
if __name__ == '__main__':
    args = parse_args()
    C = CODES3D(args.config)
    commons_db = create_engine(C.commons_db_url, echo=False, poolclass=NullPool)
    validate_args(args, commons_db)

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    logger = Logger(logfile=os.path.join(
        args.output_dir, 'codes3d.log'))
    start_time = time.time()
    tissues = parse_tissues(
        args.tissues, args.match_tissues, args.eqtl_project, commons_db)

    hic_df = parse_hic(
        args.match_tissues,
        args.include_cell_lines,
        args.exclude_cell_lines,
        args.restriction_enzymes,
        commons_db,
        args.pchic)
    
    if not args.pchic and args.match_tissues and not args.non_spatial:
        logger.write('Using Hi-C libraries to get spatial information.\t{}'.format(
            ', '.join(hic_df['library'].tolist())))
        logger.write('\neQTL tissues:\t{}\n'.format(
            ', '.join(tissues['name'].tolist())))

        upsert = input( 
            '''WARNING: We've tried to match your Hi-C and eQTL tissues above.
            Continue? [y/N]'''
        )
        if not upsert.lower() == 'y':
            print(('Use -t and -n to include specific eQTL tissues'
                   ' and Hi-C libraries'))
            sys.exit('Exiting.')


    if args.pchic and args.match_tissues and not args.non_spatial:
        logger.write('Using PCHi-C libraries to get spatial information.\t{}'.format(
            ', '.join(hic_df['library'].tolist())))
        logger.write('\neQTL tissues:\t{}\n'.format(
            ', '.join(tissues['name'].tolist())))
        
        upsert = input(
            '''WARNING: We've tried to match your PCHi-C and eQTL tissues above.
            Continue? [y/N]'''
        )
        if not upsert.lower() == 'y':
            print(('Use -t and -n to include specific eQTL tissues'
                   ' and PCHi-C libraries'))
            sys.exit('Exiting.')
    
    log_settings(args, logger)

    eqtl_project = tissues['project'].tolist()[0]
    config = configparser.ConfigParser()
    config.read(args.config)
    covariates_dir = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'COVARIATES_DIR'))
    expression_dir = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'EXPRESSION_DIR'))
    genotypes_fp = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'GENOTYPES_FP'))
    expression_table_fp = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'GENE_FP'))
    eqtl_project_db = create_engine(
        C.eqtl_db_url.format(eqtl_project.lower()), echo=False, poolclass=NullPool)
    snp_df = pd.DataFrame()
    gene_df = []
    eqtl_df = []

    if args.gene_input:
        gene_info_df = genes.get_gene_info(
            args.gene_input,
            hic_df,
            args.output_dir,
            commons_db,
            logger,
            args.pchic,
            args.suppress_intermediate_files)
        if args.non_spatial:
            map_non_spatial_eqtls(
                pd.DataFrame(),
                gene_info_df,
                tissues,
                hic_df,
                C,
                args,
                genotypes_fp,
                eqtl_project_db,
                commons_db,
                covariates_dir,
                expression_dir,
                logger,
                args.pchic)
        interactions_df = interactions.find_interactions(
            gene_info_df,
            'gencode_id',
            C.lib_dir,
            hic_df,
            args.num_processes,
            logger,
            args.pchic,
            args.suppress_intermediate_files)
        snp_df, gene_df, eqtl_df = snps.find_snps(
                interactions_df,
                gene_info_df,
                tissues,
                args.output_dir,
                C,
                genotypes_fp,
                eqtl_project_db,
                covariates_dir,
                expression_dir,
                args.pval_threshold,
                args.maf_threshold,
                args.fdr_threshold,
                args.num_processes,
                commons_db,
                logger,
                args.pchic,
                args.suppress_intermediate_files)

        if gene_df.empty or snp_df.empty or eqtl_df.empty:
            logger.write('No eQTLs found.\nProgram exiting.')
            sys.exit()
        if not args.suppress_intermediate_files:
            snp_df.to_csv(os.path.join(
                args.output_dir, 'snps.txt'), sep='\t', index=False)
            gene_df.to_csv(os.path.join(
                args.output_dir, 'genes.txt'), sep='\t', index=False)
            eqtl_df.to_csv(os.path.join(
                args.output_dir, 'eqtls.txt'), sep='\t', index=False)
    
    if args.snp_input or args.snps_within_gene:
        interactions_df = []
        gene_df = []
        gene_info_df = None
        if args.snps_within_gene:
            gene_info_df = genes.get_gene_info(
                    args.snps_within_gene,
                    hic_df,
                    args.output_dir,
                    commons_db,
                    logger,
                    args.pchic,
                    args.suppress_intermediate_files)
        snp_df = snps.get_snp(
                args.snp_input,
                gene_info_df,
                hic_df,
                args.output_dir,
                eqtl_project_db,
                C.rs_merge_arch_fp,
                logger,
                args.pchic,
                args.suppress_intermediate_files)
        if not args.suppress_intermediate_files:
            snp_df.to_csv(os.path.join(
                args.output_dir, 'snps.txt'), sep='\t', index=False)
        snp_list = snp_df['snp'].drop_duplicates().tolist()

        if args.non_spatial:
            map_non_spatial_eqtls(
                    snp_df,
                    pd.DataFrame(),
                    tissues,
                    hic_df,
                    C,
                    args,
                    genotypes_fp,
                    eqtl_project_db,
                    commons_db,
                    covariates_dir,
                    expression_dir,
                    logger,
                    args.pchic)
        else:
            gene_df, eqtl_df = map_spatial_eqtls(
                    snp_list,
                    C,
                    hic_df,
                    args,
                    logger,
                    commons_db,
                    tissues,
                    genotypes_fp,
                    eqtl_project_db,
                    covariates_dir,
                    expression_dir,
                    args.pchic
                    )
    
    eqtl_df = multi_test_correction(eqtl_df, args.multi_test.lower())
    #logger.write('  * {} eQTL associations passed FDR <= {}.'.format(
    #    len(eqtl_df[eqtl_df['adj_pval'] <= args.fdr_threshold]),
    #    args.fdr_threshold))
    logger.write('  * eQTLs mapped at MAF >= {} and pval threshold <={}.'.format(
        args.maf_threshold, args.pval_threshold))
    if not args.suppress_intermediate_files:
        eqtl_df.to_csv(os.path.join(
            args.output_dir, 'eqtls.txt'), sep='\t', index=False)
    eqtl_df = eqtl_df[eqtl_df['adj_pval'] <= args.fdr_threshold]

    if len(eqtl_df) == 0:
        logger.write('  * No significant eQTL associations found at the FDR <= {}'.format(
            args.fdr_threshold))
        logger.write('  * Program exiting.')
        sys.exit()
        
    if not args.do_not_produce_summary:
        if not args.no_afc:
            afc_start_time = time.time()
            eqtl_df = calc_afc(
                eqtl_df,
                genotypes_fp,
                expression_dir,
                covariates_dir,
                eqtl_project,
                args.output_dir,
                args.fdr_threshold,
                args.afc_bootstrap,
                args.num_processes)
        summary.produce_summary(
            args.pchic, eqtl_df, snp_df, gene_df, expression_table_fp,
            args.fdr_threshold, args.output_dir,
            args.num_processes, args.output_format, args.no_afc, logger)
    msg = 'Done.\nTotal time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60)
    logger.write(msg)
    now = datetime.datetime.now()
    now = now.strftime("%d/%m/%Y %H:%M:%S")
    logger.write(f'{now}')

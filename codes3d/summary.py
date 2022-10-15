#!/usr/bin/env python

import pandas as pd
import numpy as np
import time
import re
import os
import sys
import argparse
import configparser
import multiprocessing
from more_itertools import sliced

import codes3d


def produce_summary(
        eqtl_df, snp_df, gene_df,
        expression_table_fp, fdr_threshold, output_dir, num_processes,
        output_format, no_afc, logger, pchic=False):
    """Write final results of eQTL-eGene associations

    Args:

    Returns:
        significant_eqtls.text: A file with eQTL associations with
            adj_p_values <= FDR threshold. 
    """
    logger.write('Writing final output...')
    start_time = time.time()
    to_file = []
    global gene_exp
    gene_exp = {}
    snps_genes = []
    genes_tissues = []
    eqtl_df = eqtl_df.merge(snp_df[['snp', 'variant_id']].drop_duplicates(),
                            left_on='sid', right_on='variant_id', how='left')
    eqtl_df = eqtl_df.merge(
        gene_df[
            ['gencode_id', 'gene', 'gene_chr', 'gene_start', 'gene_end']
        ].drop_duplicates(),
        left_on='pid', right_on='gencode_id', how='left')
    eqtl_df['sid_chr'] = eqtl_df['sid_chr'].str[3:]
    gene_df['pid'] = gene_df['gencode_id'].str.split('.', expand=True)[0]
    cols = []
    if pchic:
        df = gene_df[['snp', 'gencode_id', 'cell_line','N_reads','Score']].merge(
                eqtl_df[['snp', 'gencode_id']],  on=['snp', 'gencode_id'], how='inner')
        df = df.drop_duplicates()
        grouped_df = df.groupby(['snp', 'gencode_id'])
        df['cell_lines'] = df.apply(lambda row: ', '.join(
            grouped_df.get_group((row['snp'], row['gencode_id']))[
                'cell_line'].values.tolist()), axis=1)
        df['chicago_scores'] = df.apply(lambda row: ', '.join(
            grouped_df.get_group((row['snp'], row['gencode_id']))[
                'Score'].astype(str).values.tolist()), axis=1)
        df['N_reads_per_cell_line'] = df.apply(lambda row: ', '.join(
            grouped_df.get_group((row['snp'], row['gencode_id']))[
                'N_reads'].astype(str).values.tolist()), axis=1)
        eqtl_df = eqtl_df.merge(df[
            ['snp', 'gencode_id', 'cell_lines', 'N_reads_per_cell_line', 'chicago_scores']
            ].drop_duplicates(), on=['snp', 'gencode_id'], how='left')
     
         
    else:
        logger.write("  * Computing Hi-C data...")
        df = gene_df[['snp', 'gencode_id', 'cell_line', 'interactions',
            'replicates', 'cell_line_replicates']].merge(
            eqtl_df[['snp', 'gencode_id']],  on=['snp', 'gencode_id'], how='inner')
        df['hic_score'] = df['interactions'] / df['cell_line_replicates']
        df = df.drop_duplicates()
        grouped_df = df.groupby(['snp', 'gencode_id'])
        df['cell_line_hic_scores'] = df.apply(lambda row: ', '.join(
            grouped_df.get_group((row['snp'], row['gencode_id']))[
                'hic_score'].round(2).astype(str).values.tolist()), axis=1)
        df['hic_score'] = df.apply(lambda row: grouped_df.get_group(
            (row['snp'], row['gencode_id']))[
            'hic_score'].sum().round(2), axis=1)
        df['cell_lines'] = df.apply(lambda row: ', '.join(
            grouped_df.get_group((row['snp'], row['gencode_id']))[
                'cell_line'].values.tolist()), axis=1)
        eqtl_df = eqtl_df.merge(df[
            ['snp', 'gencode_id', 'cell_lines', 'hic_score', 'cell_line_hic_scores']
            ].drop_duplicates(), on=['snp', 'gencode_id'], how='left')
    # Get SNP-gene distance
    eqtl_df['sid_chr'] = 'chr'+eqtl_df['sid_chr'].astype(str)
    eqtl_df['distance'] = eqtl_df.apply(
            lambda row: get_snp_gene_distance(row), axis=1)
    eqtl_df['interaction_type'] = eqtl_df.apply(
            lambda row: label_cis(row), axis=1)
    eqtl_df = eqtl_df.rename(
            columns={'sid_chr': 'snp_chr',
                    'sid_pos': 'snp_locus',
                    'pval': 'eqtl_pval',
                    'b': 'beta',
                    'b_se': 'beta_se'})

    if not no_afc and not pchic:
        cols = ['snp', 'gencode_id', 'gene', 'tissue', 'adj_pval',
                'log2_aFC', 'log2_aFC_lower', 'log2_aFC_upper', 'maf',
                'interaction_type', 'hic_score']
    elif no_afc and not pchic:
        cols = ['snp', 'gencode_id', 'gene', 'tissue', 'adj_pval',
                'beta', 'beta_se', 'maf', 'interaction_type', 'hic_score']

    elif not no_afc and pchic:
        cols = ['snp', 'gencode_id', 'gene', 'tissue', 'adj_pval',
                'log2_aFC', 'log2_aFC_lower', 'log2_aFC_upper', 'maf',
                'interaction_type', 'chicago_scores']
    elif no_afc and pchic:
        cols = ['snp', 'gencode_id', 'gene', 'tissue', 'adj_pval',
                'beta', 'beta_se', 'maf', 'interaction_type', 'chicago_scores']

        
    if output_format == 'short':
        eqtl_df[cols].to_csv(os.path.join(output_dir, 'significant_eqtls.txt'), sep='\t',
                             columns=cols, index=False)
        return
    eqtl_df['ref'] = eqtl_df.apply(
        lambda row: row['sid'].split('_')[2], axis=1)
    eqtl_df['alt'] = eqtl_df.apply(
        lambda row: row['sid'].split('_')[3], axis=1)
    cols += ['eqtl_pval', 'snp_chr', 'snp_locus', 'ref', 'alt',
             'gene_chr', 'gene_start', 'gene_end', 'distance']
    if output_format == 'medium':
        eqtl_df = eqtl_df[cols]
        eqtl_df.to_csv(os.path.join(output_dir, 'significant_eqtls.txt'), sep='\t',
                       index=False)
        return
    logger.write("  * Collecting gene expression rates...")
    exp_df = pd.read_table(expression_table_fp, skiprows=2, engine='c')

    # Modify columns in GTEx v8 gene expression table to match tissue names
    new_columns = []
    for column in exp_df.columns:
        column = re.sub(r"\s", '_', column)
        column = re.sub(r"['\(\)']", '', column)
        column = re.sub(r"_-_", '_', column)
        new_columns.append(column)
    exp_df.columns = new_columns

    exp_df.set_index('Name', inplace=True)

    # Get gene expression in eQTL tissue
    eqtl_df['expression'] = eqtl_df.apply(lambda row: exp_df.loc[
        row['gencode_id'], row['tissue']], axis=1)

    # Determine gene expression extremes
    extremes_df = {
        'max_expressed_tissue': exp_df.drop(['Description'], axis=1).replace(
            0, np.nan).idxmax(axis=1),
        'max_expression': exp_df.drop(['Description'], axis=1).replace(
            0, np.nan).max(axis=1),
        'min_expressed_tissue': exp_df.drop(['Description'], axis=1).replace(
            0, np.nan).idxmin(axis=1),
        'min_expression': exp_df.drop(['Description'], axis=1).replace(
            0, np.nan).min(axis=1),
    }
    extremes_df = pd.DataFrame(extremes_df)

    if pchic:
        eqtl_df = eqtl_df.merge(
            extremes_df.reset_index(), left_on=['gencode_id'], right_on=['Name'])
        cols += ['cell_lines', 'N_reads_per_cell_line', 'chicago_scores', 'expression',
                'max_expressed_tissue', 'max_expression',
                'min_expressed_tissue', 'min_expression']
    else:
        eqtl_df = eqtl_df.merge(
                extremes_df.reset_index(), left_on=['gencode_id'], right_on=['Name'])
        cols += ['cell_lines', 'cell_line_hic_scores', 'expression',
                'max_expressed_tissue', 'max_expression',
                'min_expressed_tissue', 'min_expression'] 
    eqtl_df = eqtl_df[cols]
    eqtl_df.to_csv(os.path.join(output_dir, 'significant_eqtls.txt'), sep='\t',
                   columns=cols, index=False)
    logger.write('  Time elasped: {:.2f} mins.'.format(
        (time.time() - start_time)/60))


def get_snp_gene_distance(row):
    distance = ''
    if str(row['sid_chr']) == str(row['gene_chr']):
        if row['sid_pos'] >= row['gene_start'] and \
           row['sid_pos'] <= row['gene_end']:
            distance = 0
        else:
            distance = min(abs(int(row['sid_pos']) - int(row['gene_end'])),
                           abs(int(row['gene_start']) - int(row['sid_pos'])))
    else:
        distance = 'NA'
    return distance


def label_cis(row):
    try:
        if int(row['distance']) < 1000000:
            return 'Cis'
        else:
            return 'Trans-intrachromosomal'
    except ValueError:  # Interchromosomal interactions have distance='NA'
        return 'Trans-interchromosomal'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-e', '--eqtl_files', nargs='+', default=None,
        help='''Space-separated list of one or more eQTL files
              'generated by eqtls.py. See --batches_dir''')
    parser.add_argument(
        '-b', '--batches_dir', default=None,
        help='''Filepath to the directory containing batch subdirectories
              'each containing eqtls.txt, snps.txt and genes.txt.''')
    parser.add_argument(
        '-o', '--output_dir', default='codes3d_summary',
        help='The directory in which to output results ')
    parser.add_argument(
        '--multi-test', default = 'all',
        help='''Options for BH multiple-testing: ['snp', 'tissue', 'multi'].
        'snp': corrects for genes associated with a given SNP in a given tissue.
        'tissue': corrects for all associations in a given tissue. 
        'multi': corrects for all associations across all tissues tested.''')
    parser.add_argument(
        '-f', '--fdr_threshold', type=float, default=0.05,
        help='''The FDR threshold to consider an eQTL statistically
        significant (default: 0.05).''')
    parser.add_argument(
        '--eqtl-project', type=str, default='GTEx',
        help='''The eQTL project to query. Default: GTEx.
        'use \'codes3d.py --list-eqtl-db\' to list available databases.''')
    parser.add_argument(
        "--no-afc", action="store_true", default=False,
        help='''Do not calculate allelic fold change (aFC).
        If true, eQTL beta is calculated instead. (default: False).''')
    parser.add_argument(
        '--afc-bootstrap', type=int, default=1000,
        help='Number of bootstrap for  aFC calculation (default: 1000).')
    parser.add_argument(
        '--output-format', type=str, default='long',
        help='''Determines columns of output.
        \'short\':snp|gencode_id|gene|tissue|pval|adj_pval|maf.
        \'medium\' adds:snp_chr|snp_locus|gene_chr|gene_start|
        gene_end|snp_gene_distance|interaction.
        \'long\'  adds:cell_lines|hic_scores|total_hic_score|
        gene_expression|max_expressed_tissue|max_expression|
        min_expressed_tissue|min_expression.(default: long).''')
    parser.add_argument(
        '-p', '--num_processes', type=int,
        default=min(multiprocessing.cpu_count(), 32),
        help='''The number of processes for compilation of the results
        (default: %s).''' % str(min(multiprocessing.cpu_count(), 32)))
    parser.add_argument(
        '-c', '--config',
        default=os.path.join(os.path.dirname(
            __file__), '../docs/codes3d.conf'),
        help='The configuration file to use (default: codes3d.conf)')
    parser.add_argument(
        "--suppress-intermediate-files", action="store_true", default=False,
        help='''Do not produce intermediate files. These can be used to run the 
        pipeline from an intermediate stage in the event of interruption
        (default: False)''')
    args = parser.parse_args()
    if (not args.eqtl_files) and (not args.batches_dir):
        sys.exit('Required --eqtl-files or --batches_dir missing.')
    start_time = time.time()
    config = configparser.ConfigParser()
    config.read(args.config)
    eqtl_project = args.eqtl_project
    covariates_dir = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'COVARIATES_DIR'))
    expression_dir = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'EXPRESSION_DIR'))
    genotypes_fp = os.path.join(
        os.path.dirname(__file__), config.get(eqtl_project, 'GENOTYPES_FP'))
    expression_table_fp = os.path.join(
        os.path.dirname(__file__), config.get(args.eqtl_project, 'GENE_FP'))
    snp_files = []
    gene_files = []
    eqtl_files = []
    if args.eqtl_files:
        snp_files = [os.path.join(os.path.dirname(fp), 'snps.txt')
                     for fp in args.eqtl_files]
        gene_files = [os.path.join(os.path.dirname(fp), 'genes.txt')
                      for fp in args.eqtl_files]
        eqtl_files = args.eqtl_files
    elif args.batches_dir:
        snp_files = [os.path.join(args.batches_dir, batch_dir,  'snps.txt')
                     for batch_dir in os.listdir(args.batches_dir) if os.path.exists(
                             os.path.join(args.batches_dir, batch_dir, 'snps.txt'))]
        gene_files = [os.path.join(args.batches_dir, batch_dir,  'genes.txt')
                      for batch_dir in os.listdir(args.batches_dir) if os.path.exists(
                             os.path.join(args.batches_dir, batch_dir, 'genes.txt'))]
        eqtl_files = [os.path.join(args.batches_dir, batch_dir,  'eqtls.txt')
                      for batch_dir in os.listdir(args.batches_dir) if os.path.exists(
                             os.path.join(args.batches_dir, batch_dir, 'eqtls.txt'))]
    snp_df = codes3d.parse_intermediate_files(
        snp_files, args.output_dir, 'snps')
    gene_df = codes3d.parse_intermediate_files(
        gene_files, args.output_dir, 'genes')
    eqtl_df = codes3d.parse_intermediate_files(
        eqtl_files, args.output_dir, 'eqtls', args.multi_test)
    if eqtl_df.empty or snp_df.empty or gene_df.empty:
        sys.exit('Missing eqtls.txt, or snps.txt or genes.txt files.')
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    c = codes3d.CODES3D(args.config)
    logger = codes3d.Logger(logfile=os.path.join(
        args.output_dir, 'codes3d.log'))
    if not args.suppress_intermediate_files:
        logger.write('Writing intermediary files...')
        snp_df.to_csv(os.path.join(
            args.output_dir, 'snps.txt'), sep='\t', index=False)
        gene_df.to_csv(os.path.join(
            args.output_dir, 'genes.txt'), sep='\t', index=False)
        eqtl_df.to_csv(os.path.join(
            args.output_dir, 'eqtls.txt'), sep='\t', index=False)
    eqtl_df = eqtl_df[eqtl_df['adj_pval'] <= args.fdr_threshold].sort_values(by=['sid']).reset_index(drop=True)
    if not args.no_afc:
        chunk_size = 1000
        index_slices = sliced(range(len(eqtl_df)), chunk_size)
        afc_df = []
        for index_slice in index_slices:
            afc_df.append(
                codes3d.calc_afc(
                    eqtl_df.iloc[index_slice], genotypes_fp, expression_dir, covariates_dir,
                    eqtl_project, args.output_dir, args.fdr_threshold, args.afc_bootstrap,
                    args.num_processes)
                )
        eqtl_df = pd.concat(afc_df)
    produce_summary(
        eqtl_df, snp_df, gene_df, expression_table_fp,
        args.fdr_threshold, args.output_dir,
        args.num_processes, args.output_format, args.no_afc, logger)
    logger.write('Time elapsed: {:.2f}'.format((time.time()-start_time)/60))

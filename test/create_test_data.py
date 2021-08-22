import pandas as pd
import os
from sqlalchemy import create_engine
import codes3d
import snps
import sys


def get_snp_list(test_dir):
    '''
    rs1794265 # biallelic
    rs1689510 # at start of MboI fragment; ensures fragments are 0-based
    '''
    return pd.read_csv(
        os.path.join(test_dir, 'test_snps.txt'),
        names=['snp']
    )


def get_merged_snp(test_dir):
    '''
    rs116755193 # has been merged into rs34660260 
    '''
    return pd.read_csv(
        os.path.join(test_dir, 'merged_snps.txt'),
        names=['snp']
    )


def get_snp_rs_df(test_dir):
    df = pd.read_csv(
        os.path.join(test_dir, 'snps.txt'),
        sep='\t')
    df = df.rename(columns={'snp': 'rsid'})
    return df[
        ['id', 'rsid', 'chr', 'locus', 'variant_id']
    ].drop_duplicates().reset_index(drop=True).sort_values(
        by=['id'], ascending=True)


def get_snp_fragments(test_dir):
    df = pd.read_csv(
        os.path.join(test_dir, 'snps.txt'),
        sep='\t')
    df = df.rename(
        columns={
            'snp': 'rsid',
            'fragment': 'frag_id'})
    return df.sort_values(
        by=['id', 'chr', 'frag_id', 'enzyme']
    ).reset_index(drop=True)

def create_output_dir(test_dir):
    output_dir = os.path.join(test_dir, 'output_dir')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    return output_dir

test_dir = os.path.join(os.path.dirname(__file__),
                        "../test/")
config_fp = os.path.join(os.path.dirname(__file__),
                         "../docs/codes3d.conf")
c = codes3d.CODES3D(config_fp)
db = create_engine(c.postgres_commons, echo=False)
enzymes = ['MboI', 'HindIII', 'DpnII']

output_dir = create_output_dir(test_dir)
snp_list = get_snp_list(test_dir)
snp_rsid_df = get_snp_rs_df(test_dir)
snp_fragment_df = get_snp_fragments(test_dir)
merged_snps_fp = os.path.join(test_dir, 'merged_snps.txt') 

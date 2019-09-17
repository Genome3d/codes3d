#!/usr/bin/env python

import os
import sqlite3
import pandas as pd
import progressbar
import argparse
import configparser
import sys
import shutil
import codes3d
"""
At first sqlite3 threw a 'database is locked' exception. Solved this by adding
'nobrl' to options in mount.sh. This in turn led to a 'FileNotFoundError' when
importing pandas. Using `cd $(pwd)` command resolved this.
"""


def process(gtex_zip_fp, done_dir):
    print('Creating index for ...')
    for tissue_file in os.listdir(gtex_zip_fp):
        if tissue_file.endswith('allpairs.txt.gz'):
            tissue = tissue_file[:len(tissue_file)-16]
            print('\n\t ', tissue)
            conn = sqlite3.connect(os.path.join(gtex_zip_fp, tissue+'.db'))
            conn.text_factory = str
            cur = conn.cursor()
            cur.execute("""CREATE TABLE associations (
            gene_id TEXT, 
            variant_id TEXT, 
            tss_distance INTEGER, 
            ma_samples INTEGER, 
            ma_count INTEGER, 
            maf REAL, 
            pval_nominal REAL, 
            slope REAL, 
            slope_se REAL)""")
            cur.execute(
                "CREATE INDEX IF NOT EXISTS id ON associations (gene_id,variant_id)")
            df = pd.read_csv(os.path.join(gtex_zip_fp, tissue_file),
                             compression='gzip',
                             delimiter='\t')
            bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            for index, association in df.iterrows():
                cur.execute("INSERT INTO associations VALUES(?,?,?,?,?,?,?,?,?);",
                            [association['gene_id'],
                             association['variant_id'],
                             association['tss_distance'],
                             association['ma_samples'],
                             association['ma_count'],
                             association['maf'],
                             association['pval_nominal'],
                             association['slope'],
                             association['slope_se']])
                bar.update(index)
            conn.commit()
            cur.close()
            conn.close()

            #shutil.move(os.path.join(gtex_zip_fp, tissue_file), done_dir)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', required=True,
                        help='Directory containing all GTEx analysis')
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file specifying the location \
                        of the CoDeS3D library (default: docs/codes3d.conf).")
    parser.add_argument('-o', '--output_dir', default='lib/',
                        help='Directory to put tissue database files.')
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    eqtl_data_dir = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "EQTL_DATA_DIR"))

    if not os.path.isdir(eqtl_data_dir):
        os.mkdir(eqtl_data_dir)
    codes3d.build_local_eqtl_index(args.input, eqtl_data_dir)

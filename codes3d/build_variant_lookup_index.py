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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a variant lookup table from GTEx. \
    Used in local GTEx query to convert SNP IDs to rsIDs.')
    parser.add_argument('-i', '--input', required=True,
                        help='Variant lookup table from GTEx')
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
    snp_dict_fp = os.path.join(os.path.dirname(__file__),
                               config.get("Defaults", "SNP_DICT_FP"))

    codes3d.build_variant_lookup_index(args.input, snp_dict_fp)

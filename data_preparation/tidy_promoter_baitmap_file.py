#!/usr/bin/env python 

import pandas as pd
import argparse
import time
import os
import sys


def tidy_promoter_file(baits,out_file):
    
    cols = ['chr', 'start', 'end', 'frag_id', 'genes']
    bait_df = pd.read_csv(baits, sep="\t", names = cols).drop_duplicates()
    bait_df['genes'] = bait_df['genes'].str.replace(', ', ',')
    bait_df['genes'] = bait_df['genes'].str.split(',')
    #print(bait_df)
    bait_df_new = bait_df.explode('genes')
    bait_df_new = bait_df_new.dropna()
    bait_df_new = bait_df_new.drop_duplicates()
    bait_df_new.sort_values(by = 'genes', inplace=True, ascending=True)
    print(bait_df_new)

    bait_df_new.to_csv(out_file, sep='\t', header=True, index=False)

def parse_args():

    parser = argparse.ArgumentParser(
            description="Tidy up promoter file for CoDeS3D PRO")

    parser.add_argument(
            '-f', '--bait-file', required = True,
            help = '''Pass a promoter information file with an extension .baitmap 
            which was used to run CHiCAGO experiment''')
    
    parser.add_argument(
            '-o', '--out-fp', required = True,
            help = '''Provide absolute path to an output file''')

    return parser.parse_args()


if __name__ == '__main__':

    args =  parse_args()
    start_time = time.time()

    if not (args.bait_file or args.out_fp):
        message = '''One or more of the required parameters are missing.
        Use tidy_promoter_baitmap_file.py -h for more details.'''
        print(message)
        sys.exit()

    tidy_promoter_file(args.bait_file, args.out_fp)



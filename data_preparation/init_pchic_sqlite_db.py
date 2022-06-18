#!/usr/bin/env python

import pandas as pd
import sqlite3
from sqlite3 import Error
import os
import sys
import argparse
import time

def read_digested_genome(genome):

    cols = ['fchr','fstart','fend','fid']
    if genome.endswith('.rmap'):
        genome_df = pd.read_csv(genome, sep = "\t", header = None, names=cols)
        #print(genome_df)
    return genome_df


def read_chicago_output_ibeds(inputs):
    
    inter_cols = ['b_chr','b_start','b_end','b_name','oe_chr','oe_start','oe_end','oe_name','n_reads','score']
    input_df = pd.read_csv(inputs, sep = "\t", skiprows=1, names=inter_cols)
    return input_df


def build_pchic_db(inputs, genome, tissue_dirs, reps_dirs, out_dir, db_tblname):
    
    genome_df = read_digested_genome(genome)
    input_df = read_chicago_output_ibeds(inputs)

    bait_df = pd.merge(input_df, genome_df, how='inner', left_on=['b_chr','b_start', 'b_end'],
                            right_on=['fchr','fstart', 'fend'], sort = False)
    bait_short_df = bait_df[['b_chr','b_start','b_end','b_name','fid']]
    bait_short_df.columns = ['p_chr','p_start','p_end','p_name','p_fid']
 
    oe_df = pd.merge(input_df, genome_df, how='inner', left_on=['oe_chr','oe_start', 'oe_end'],
                          right_on=['fchr','fstart', 'fend'], sort = False)
    oe_short_df = oe_df[['oe_chr','oe_start','oe_end','oe_name','fid','n_reads','score']]
    oe_short_df.columns = ['oe_chr','oe_start','oe_end','oe_name','oe_fid','n_reads','score']

    bait_oe_df = pd.concat([bait_short_df, oe_short_df], axis=1, sort=False)

    tissue_dir_basename = os.path.splitext(os.path.basename(tissue_dirs))[0]
    rep_dir_basename = os.path.splitext(os.path.basename(reps_dirs))[0]
    db_file = rep_dir_basename+'.db'
    out_db_fp = os.path.join(out_dir,tissue_dir_basename)
    
    if not os.path.exists(out_db_fp):
        os.makedirs(out_db_fp, exist_ok=True)
        os.chmod(out_db_fp, 0o777)
    
    db_fp = os.path.join(out_db_fp,db_file)
    conn = None
    with sqlite3.connect(db_fp) as conn:
        cursor = conn.cursor()
        print("Creating database {}".format(db_file))
        bait_oe_df.to_sql(db_tblname, if_exists='replace', con=conn, index=False)
        print("Creating index on tables in database {}...".format(db_file))
        cursor.execute('''CREATE INDEX idx_{}_baits
        ON {} (p_fid)'''.format(rep_dir_basename, db_tblname))
        cursor.execute('''CREATE INDEX idx_{}_oes
        ON {} (oe_fid)'''.format(rep_dir_basename, db_tblname))
        cursor.execute('''CREATE INDEX idx_{}_pfid_oefid 
        ON {} (p_fid,oe_fid)'''.format(rep_dir_basename, db_tblname))
    print("Done!")

 
def parse_args():
    
    parser = argparse.ArgumentParser(
             description = "Build SQLite libraries of PCHi-C data for CoDeS3D.")
    
    parser.add_argument(
           '-r', '--rmap-file', required = True,
           help = '''Bowtie digested hg38 build genome file with an extension (.rmap)
                  which was used to run CHiCAGO experiment.''')

    parser.add_argument(
           '-d', '--pchic-dir', required = True,
           help = '''Provide a superdirectory (i.e. parent) containing tissues subdirectories 
                  that contain replicates subdirectories (e.g. parent/tissues/replicates_{1..n}/) 
                  where chicago output .ibed file is residing.''')

    parser.add_argument(
           '-o', '--output_dir', required = True,
           help = '''Provide path to the output directory to save PCHi-C databases.
                  The script will automatically create an output directory if the provided
                  directory doesn't exist in the path''')

    return parser.parse_args()

if __name__ == '__main__':
    
    args = parse_args()
    start_time = time.time()
    
    if not (args.rmap_file or args.pchic_dir) and args.output_dir:
        message = '''Missing --rmap-file, --pchic-dir and --output_dir. 
                  See init_pchic_sqlite_db.py -h for more details.'''
        print(message)
        sys.exit()

    if args.rmap_file == None:    
        message = '''Missing parameter... Genome fragment file with extension .rmap is required.
                   For more details init_pchic_sqlite_db.py -h'''
        print(message)
        sys.exit()
    else:
        print("\nUsing rmap file in: '%s'\n" %args.rmap_file)

    if args.pchic_dir == None: 
        message = '''Missing parameter... Directory containing .ibed files from Chicago is required.
                  For more details init_pchic_sqlite_db.py -h'''
        print(message)
        sys.exit()

    if args.output_dir == None:

        message = '''Missing parameter... Output directory to save databases is required.
                 For more details init_pchic_sqlite_db.py -h'''
        print(message)
        sys.exit()    
    
    if os.path.exists(args.output_dir): #check is directory exists
        out_dir = args.output_dir
        print("PCHi-C databases will be stored in the output directory: '%s'\n" %out_dir)
    else:
        out_dir = os.makedirs(args.output_dir, exist_ok = True)
        print("PCHi-C databases will be stored in the output directory: '%s'\n" %out_dir)
        
    db_tblname = 'interactions'
    for tissues in os.listdir(args.pchic_dir):
        #print(tissues)
        tissue_dirs = os.path.join(args.pchic_dir, tissues)
        #print(tissue_dirs)
        for reps in os.listdir(tissue_dirs):
            #print(reps)
            reps_dirs = os.path.join(tissue_dirs, reps) 
            #print(reps_dirs)
            for files in os.listdir(reps_dirs):
                files_path = os.path.join(reps_dirs, files)
                #print(files_path)
                if not files_path.endswith('.ibed'):
                    pass
                else:
                    input_fp = os.path.join(reps_dirs,files_path)
                    build_pchic_db(input_fp, 
                                   args.rmap_file,
                                   tissue_dirs, 
                                   reps_dirs, 
                                   out_dir,
                                   db_tblname)
    print("All done.")
    print('Total time elapsed: {:.2f} mins'.format((time.time()-start_time)/60))

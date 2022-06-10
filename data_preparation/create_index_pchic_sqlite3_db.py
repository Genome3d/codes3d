#!/usr/bin/env python

import pandas as pd
import sqlite3
import sys
import os


def create_index_db(db_path):
    
    for path in db_path:
        conn = None
        with sqlite3.connect(path) as conn:
            cursor = conn.cursor()
            print("Creating index on tables in database...")
            cursor.execute('''CREATE INDEX idx_{}_{} ON 
                interactions (p_fid,oe_fid)'''.format('pfid','oefid'))

    return 

db_fp = '/mnt/projects/codes3d/lib/pchic/HindIII/pchic_data/'

db_path = []
for cell_lines in os.listdir(db_fp):
    cell_lines_fp = os.path.join(db_fp,cell_lines)
    #print(cell_lines_fp)
    for dbs in os.listdir(cell_lines_fp):
        dbs_fp = os.path.join(cell_lines_fp,dbs)
        db_path.append(dbs_fp)

create_index_db(db_path)

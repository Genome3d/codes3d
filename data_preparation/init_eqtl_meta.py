#! /usr/bin/env python

import os
import pandas as pd
import sys
from sqlalchemy import create_engine

"""
Create PostgreSQL table of list of eQTL projects. 
Requires a meta_eqtls.txt
"""
eqtl_fp = os.path.join(os.path.dirname(__file__), '../lib/meta_eqtls.txt')
eqtls = pd.read_csv(eqtl_fp, sep='\t')
db = create_engine(
    'postgres://codes3d:codes3d@127.0.0.1/codes3d_commons', echo=False)
eqtls.to_sql('meta_eqtls', con=db, if_exists='replace', index=False)

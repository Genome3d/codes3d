#! /usr/bin/env python

import os
import pandas as pd
import sys
from sqlalchemy import create_engine

if __name__ == '__main__':
    '''
    Create table for meta info on PCHi-C libraries.
    To be run after placing PCHi-C libraries in
    '../lib/pchic/<enzyme>/pchic_data/<tissue or cells>/<replicate>.db'
    '''
    pchic_fp = os.path.join(os.path.dirname(__file__), '/mnt/projects/codes3d/lib/pchic')
    enzyme_dict = []
    for enzyme in os.listdir(pchic_fp):
        enzyme_fp = os.path.join(pchic_fp, enzyme)
        for library in os.listdir(os.path.join(enzyme_fp, 'pchic_data')):
            library_fp = os.path.join(enzyme_fp, 'pchic_data', library)
            #for library in os.listdir(os.path.join(project_fp)):
                #library_fp = os.path.join(project_fp,library)
            for rep in os.listdir(library_fp):
                enzyme_dict.append((rep.split('.')[0], enzyme, library))
    replicates = pd.DataFrame(enzyme_dict, columns=['replicate', 'enzyme', 'library'])
    libraries = replicates.groupby(
            ['library','enzyme']).size().reset_index(name='rep_count')
    from_file = pd.read_csv('/mnt/projects/users/sgok603/codes3d_sreemol/lib/meta_info/meta_pchic.txt', sep="\t")
    libraries = libraries.merge(from_file, how='left', on=['library','enzyme']).drop_duplicates()
    db = create_engine(
            'postgresql://codes3d:codes3d@127.0.0.1/codes3d_commons', echo=False)
    libraries.drop_duplicates().to_sql('meta_pchic', con=db, if_exists='replace')

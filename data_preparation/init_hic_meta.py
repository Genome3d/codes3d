#! /usr/bin/env python

import os
import pandas as pd
import sys
from sqlalchemy import create_engine


if __name__=='__main__':
    """Create table for meta info on Hi-C libraries.
    To be run after placing Hi-C libraries in 
    '../lib/hic/<Enzyme>/hic_data/<Cell line>/<Replicate>.db'
    """
    hic_fp = os.path.join(os.path.dirname(__file__), '../lib/hic')
    enzyme_dict = []
    for enzyme in os.listdir(hic_fp):
        enzyme_fp = os.path.join(hic_fp, enzyme)
        for library in os.listdir(os.path.join(enzyme_fp, 'hic_data')):
            library_fp = os.path.join(enzyme_fp, 'hic_data', library)
            for rep in os.listdir(library_fp):
                enzyme_dict.append((rep.split('.')[0], enzyme, library))
    replicates = pd.DataFrame(
        enzyme_dict, columns=['replicate', 'enzyme', 'library'])
    libraries = replicates.groupby(
        ['library', 'enzyme']).size().reset_index(name='rep_count')
    from_file = pd.read_csv('/mnt/projects/codes3d/lib/meta_info/meta_hic.txt', sep='\t')
    libraries = libraries.merge(
        from_file, how='left', on=['library', 'enzyme']
    ).drop_duplicates()

    db = create_engine(
        'postgres://codes3d:codes3d@127.0.0.1/codes3d_commons', echo=False)
    libraries.drop_duplicates().to_sql('meta_hic', con=db, if_exists='replace')

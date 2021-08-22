#! /usr/bin/env python
import argparse
import sys
import os
import pandas as pd
#import subprocess


def create_sample_participant_lookup(data_dir):
    sampleDS_fp = os.path.join(data_dir,
                               'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
    subjectDS_fp = os.path.join(data_dir,
                                'GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
    if os.path.isfile(sampleDS_fp):
        sampleDS = pd.read_csv(sampleDS_fp, sep='\t', index_col='SAMPID')
    else:
        print('Could not find {} in the data_dir you entered'
              .format('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'))
        sys.exit()
    if os.path.isfile(subjectDS_fp):
        subjectDS = pd.read_csv(subjectDS_fp, sep='\t', index_col='SUBJID')
    else:
        print('Could not find {} in the data_dir you entered'
              .format('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt'))
        sys.exit()
    sample_participant = []
    for subj_id, subject in subjectDS.iterrows():
        for samp_id, sample in sampleDS.iterrows():
            if subj_id in samp_id:
                sample_participant.append([samp_id, subj_id])
                print(subj_id, samp_id)
    df = pd.DataFrame(sample_participant,
                      columns=['sample_id', 'participant_id'])
    df.to_csv(os.path.join(args.data_dir, 'sample_participant_lookup.txt'),
              sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate normalized expression BED files for eQTL analyses')
    # parser.add_argument(
    #    'tpm_gct', help='GCT file with expression in normalized units, e.g., TPM or FPKM')
    #parser.add_argument('counts_gct', help='GCT file with read counts')
    #parser.add_argument('annotation_gtf', help='GTF annotation')
    parser.add_argument('data_dir',
                        help='Directory containing GTEx v8 files')
    args = parser.parse_args()
    if not os.path.isfile(os.path.join(args.data_dir, 'sample_participant_lookup.txt')):
        print('Creating sample_participant_lookup file...')
        create_sample_participant_lookup(args.data_dir)
    else:
        print('data_dir has sample_participant_lookup file')

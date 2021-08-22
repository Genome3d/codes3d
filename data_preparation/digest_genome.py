#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO
from Bio.Restriction import Restriction_Dictionary, RestrictionBatch
import csv
import sqlite3
import sys
import pandas as pd


def digest_genome(
        genome_fp,
        restriction_enzyme,
        output_dir,
        linear=False):
    base_fp = os.path.basename(genome_fp)
    if '.' in base_fp:
        base_fp = '{}.{}.fragments.bed'.format(
            base_fp[:base_fp.rfind('.')], restriction_enzyme)
    else:
        base_fp = '{}.{}.fragments.bed'.format(base_fp, restriction_enzyme)
    base_fp = os.path.join(output_dir, base_fp)
    if os.path.isfile(base_fp):
        overwrite = input(
            'WARNING: Overwriting existing fragment BED {}. Continue? [y/N]'
            .format(base_fp))
        if not overwrite.lower() == 'y':
            print("Did not overwrite existing fragment BED.")
            return
        os.remove(base_fp)
    print("Digesting")
    genome = None
    if "fasta" in genome_fp or "fa" in genome_fp:
        genome = SeqIO.parse(open(genome_fp, "rU"), format='fasta')
    else:
        genome = SeqIO.parse(open(genome_fp, "rU"), format='genbank')
    for chromosome in genome:
        print('{}\t{}'.format(chromosome.id,  len(chromosome.seq)))
        # Digest the sequence data and return the cut points
        enzyme = RestrictionBatch([restriction_enzyme])
        for enzyme, cutpoints in enzyme.search(
                chromosome.seq, linear=linear).items():
            if len(cutpoints) == 0:
                print('No restriction sites found for {}'.format(
                    chromosome.id))
                continue
            df = pd.DataFrame(cutpoints, columns=['cutpoint'])
            df['end'] = df.cutpoint-1
            df['start'] = df.end - (df.cutpoint.diff())
            df.loc[0, 'start'] = 0
            df['start'] = df['start'].astype('Int64')
            if len(df) > 1:
                last_fragment = pd.DataFrame(
                    {'start': [df.loc[len(df)-1, 'end']],
                     'end': [len(chromosome.seq)],
                     'cutpoint': [-1]})
                df = df.append(last_fragment, ignore_index=True)
            else:
                df.loc[len(df)-1, 'end'] = len(chromosome.seq)
            df['frag_id'] = df.index
            # chromosome has 'chr'
            accession = chromosome.id
            version = ''
            if "." in chromosome.id:
                accession, version = chromosome.id.split(".")
            if not accession.startswith("chr"):
                accession = "chr" + accession
            df['chr'] = accession
            df[['chr', 'start', 'end', 'frag_id']].to_csv(
                base_fp, index=False, sep='\t',
                mode='a', header=None)


def build_fragment_index(fragment_fp, output_db):
    if not output_db:
        if not fragment_fp.rfind('.') == -1:
            output_db = fragment_fp[:fragment_fp.rfind('.')] + ".db"
        else:
            output_db = fragment_fp + ".db"
    if os.path.isfile(output_db):
        overwrite = input(
            "WARNING: Overwriting existing fragment database %s. Continue? [y/N] " %
            output_db)
        if not overwrite.lower() == 'y':
            print("Did not overwrite existing fragment database.")
            return
        os.remove(output_db)
    if not os.path.isdir(os.path.dirname(output_db)):
        os.path.makedirs(os.path.dirname(output_db))
    fragment_index_db = sqlite3.connect(output_db)
    fragment_index = fragment_index_db.cursor()
    fragment_index.execute(
        "CREATE TABLE fragments (chr text, start integer, end integer, fragment integer)")
    fragment_index.execute("CREATE INDEX f_index ON fragments (chr,fragment)")

    with open(fragment_fp, 'r') as fragments_bed:
        for line in fragments_bed:
            fragment = line.strip().split('\t')
            fragment_index.execute("INSERT INTO fragments VALUES (?,?,?,?)", [
                                   fragment[0][fragment[0].find("chr") + 3:], int(fragment[1]), int(fragment[2]), fragment[3]])
    fragment_index_db.commit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            'Digest a genome with a restriction enzyme. ',
            'Creates a BED file where the fourth column denotes the restriction',
            ' fragment number of the fragment specified'))
    parser.add_argument(
        "-g", "--genome",
        help="Genome in genbank or fasta format.")
    parser.add_argument(
        "-l", "--linear", action="store_true",
        help="Is the genome linear? (default: False)")
    parser.add_argument(
        "-e", "--enzyme", nargs='+',
        help=("The restriction enzyme with which to fragment the genome, ",
              "for example, `-e MspI`"))
    parser.add_argument(
        "-o", "--output_dir",
        help=("Output directory for fragment BED and database files. ",
              "(default: the directory of input genome.)"))
    parser.add_argument(
        "-b", "--do_not_index", action="store_true",
        help=("Suppress building of fragment index from resultant ",
              "fragmented genome."))
    parser.add_argument(
        "-d", "--output_db",
        help=("Output file path for fragment BED file and database ",
              "(default: the name of the input file, with the extension ",
              "\".db\"/\".bed\", respectively.)"))
    parser.add_argument(
        "-z", "--list_enzymes",
        help="List the available enzymes for digestion.", action="store_true")
    args = parser.parse_args()

    enzyme_ref = Restriction_Dictionary.rest_dict
    if args.list_enzymes:
        for enzyme in enzyme_ref:
            print('{}\t{}'.format(enzyme, enzyme_ref[enzyme]['site']))
        sys.exit()
    if not args.genome:
        sys.exit('genome file is required.')
    enzymes = []
    if args.enzyme:
        enzymes_lower = [e.lower() for e in list(enzyme_ref.keys())]
        enzymes = [list(enzyme_ref.keys())[enzymes_lower.index(
            e.lower())] for e in args.enzyme if e.lower() in enzymes_lower]
    output_dir = None
    if not args.output_dir:
        output_dir = os.path.dirname(args.genome)
        if output_dir == '':
            output_dir = './'
    else:
        output_dir = args.output_dir
    if os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    for enzyme in enzymes:
        digest_genome(
            args.genome, enzyme, output_dir,
            args.linear)

        # if not args.do_not_index:
        #    build_fragment_index(base_fp, output_db)

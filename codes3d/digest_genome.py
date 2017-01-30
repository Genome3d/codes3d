#!/usr/bin/python

import os
import sys
import csv
import codes3d

from argparse import ArgumentParser
from Bio import SeqIO
from Bio import Restriction
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import Restriction_Dictionary

def digest(args):
    print "Digesting"
    enzyme = RestrictionBatch([args.enzyme])
    if "fasta" in args.genome or "fa" in args.genome:
        genome = SeqIO.parse(open(args.genome,"rU"), format='fasta')
    else:
        genome = SeqIO.parse(open(args.genome,"rU"), format='genbank')

    with csv.writer(open(args.output,"w"),delimiter=",") as output:
        output.writerow(['accession','version','start','end','sequence_start','sequence_end'])

        for chromosome in genome:
            print chromosome.id, chromosome.name
            #Digest the sequence data and return the cut points
            for enzyme, cutpoints in enzyme.search(chromosome.seq, linear=args.linear).items():
                cutpoints.insert(0,0)
                #Covert cut points to fragments
                for index, point in enumerate(cutpoints):
                    #Adjust for start and end of sequence and the offset for the cutpoint
                    #ATTENTION: This will only work for MspI I will have to spend some time to make it compatible with 
                    #           any restriction enzyme such as ones with blunt ends
                    #startpoint = 0 if cutpoints[index] - 1 < 0 else cutpoints[index] - 1
                    startpoint = 0 if cutpoints[index] - 1 < 0 else cutpoints[index]
                    endpoint = len(chromosome.seq) if index+1 >= len(cutpoints) else cutpoints[index+1] - 1

                    accession = chromosome.id
                    version = ''
                    if "." in chromosome.id:
                        accession, version = chromosome.id.split(".")
                    
                    output.writerow([accession, version, startpoint, endpoint, chromosome.seq[startpoint:startpoint+75],chromosome.seq[endpoint-75:endpoint]])
    if args.build_fragment_index:
        if not args.output_bed_db:
            if not args.output.rfind('.') == -1:
                args.output_bed_db = args.output[:args.output.rfind('.')]
            else:
                args.output_bed_db = args.output
        codes3d.build_fragment_index(args.output,args.output_bed_db)
                       
if __name__ == '__main__':
    parser = ArgumentParser(description="Digest a genome with a restriction enzyme.")
    parser.add_argument("-g","--genome", help = "Genome in genbank or fasta format.")
    parser.add_argument("-l","--linear", action="store_true", help = "Is the genome linear? (default: False)", default=False)
    parser.add_argument("-e","--enzyme", help = "The restriction enzyme with which to fragment the genome, for example, `-e MspI`")
    parser.add_argument("-o","--output", help = "Output file path.")
    parser.add_argument("-b","--build_fragment_index", action="store_true", default=False, help = "Build fragment index and BED file from resultant fragmented genome.")
    parser.add_argument("-p","--output_bed_db", help = "Output file path for fragment BED file and database (default: the name of the input file, with the extension \".db\"/\".bed\", respectively.)")

    parser.add_argument("-z","--list_enzymes", help= "List the available enzymes for digestion.", action="store_true")
    args = parser.parse_args()

    if args.list_enzymes:
        for enzyme in Restriction_Dictionary.rest_dict.keys():
            print enzyme, Restriction_Dictionary.rest_dict[enzyme]['site']
    else:
        digest = Digest(args)
        digest.start()


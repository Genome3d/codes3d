#!/usr/bin/env python

import argparse,codes3d

from Bio.Restriction import Restriction_Dictionary

                   
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Digest a genome with a restriction enzyme. Creates a BED file where the fourth column denotes the restriction fragment number of the fragment specified")
    parser.add_argument("-g","--genome", help = "Genome in genbank or fasta format.")
    parser.add_argument("-l","--linear", action="store_true", help = "Is the genome linear? (default: False)")
    parser.add_argument("-e","--enzyme", help = "The restriction enzyme with which to fragment the genome, for example, `-e MspI`")
    parser.add_argument("-o","--output", help = "Output file path.")
    parser.add_argument("-b","--do_not_index", action="store_true", help = "Suppress building of fragment index from resultant fragmented genome.")
    parser.add_argument("-d","--output_db", help = "Output file path for fragment BED file and database (default: the name of the input file, with the extension \".db\"/\".bed\", respectively.)")
    parser.add_argument("-z","--list_enzymes", help= "List the available enzymes for digestion.", action="store_true")
    args = parser.parse_args()

    if args.list_enzymes:
        for enzyme in Restriction_Dictionary.rest_dict.keys():
            print(enzyme, Restriction_Dictionary.rest_dict[enzyme]['site'])
    else:
        codes3d.digest_genome(args.genome,args.enzyme,args.output,args.output_db,args.do_not_index,args.linear)


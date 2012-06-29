#!/usr/bin/python

import os
import sys
import math

from snptools import SNPLine
from check_snps import read_snp_data

###
def main(argv):

    snp_data = read_snp_data(open(argv[0]))
    outfile = 'ped-0-snps.iit.txt'
    if len(argv) >= 2:
        outfile = argv[1]

    fo = open(outfile, 'w')

    for snpline in sorted(snp_data.items()):
        snpline = snpline[1]
        snp = (snpline.lineid, snpline.seqid, snpline.pos,
               snpline.refbase, snpline.mutation)
        fo.write('>ped%i %s:%i %c%c\n' % snp)

    fo.close()


    return None

if __name__ == '__main__': main(sys.argv[1:])

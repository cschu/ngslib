#!/usr/bin/python

import os
import sys
import math

from snptools import SNPLine
from check_snps import read_snp_data

###
def main(argv):

    snp_data = read_snp_data(open(argv[0]))

    fo = open('ped-0-snps.iit.txt', 'w')

    for snpline in sorted(snp_data.items()):
        snpline = snpline[1]
        snp = (snpline.lineid, snpline.seqid, snpline.pos,
               snpline.refbase, snpline.mutation)
        fo.write('>ped%i %s:%i %c%c\n' % snp)

    fo.close()


    return None

if __name__ == '__main__': main(sys.argv[1:])

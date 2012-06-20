#!/usr/bin/python

import os
import re
import sys
import math

import CodingDictionary as cd
from snptools import SNPLine



###
def main(argv):

    fi = open(argv[0])
    coding_d = cd.CodingDictionary(open(argv[1]),
                                   {'exon': 1, 'CDS': 2})
    
    # Test block
    # print 13538 in coding_d['1']
    # print coding_d['1'][13538]
    # print coding_d.query_position('1', 13538, 1)
    
    gffout = open('ped-0_snps_in_exons.gff', 'w')
    tabout = open('ped-0_snps_in_exons.txt', 'w')
    i = 1
    while True:
        line = fi.readline()
        if not line: break
        snpline = SNPLine(line, lineid=i)
        if coding_d.query_position(snpline.seqid[-1],
                                   snpline.pos,
                                   1):
            gffout.write('%s\n' % snpline.get_gffstring())
            tabout.write('%s\n' % str(snpline))
            i += 1
        pass
    gffout.close()
    tabout.close()

    fi.close()


    return None

if __name__ == '__main__': main(sys.argv[1:])

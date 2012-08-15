#!/usr/bin/env python

'''
Created on Aug 15, 2012

@author: Chris
'''
import sys

import pysam

from snptools import SNPDict
import gff_helpers

MIN_NREADS = 3


# region = (gffline[0], int(gffline[3]) - 1, 
#                  int(gffline[4]) - 1, gffline[6], gene_id)

# self[(snpline.seqid, snpline.pos)] = snpline Chr%

"""def binary_search(a, x, lo=0, hi=None):
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        midval = a[mid]
        if midval < x:
            lo = mid+1
        elif midval > x: 
            hi = mid
        else:
            return mid
    return -1

"""
def find_region(pos, regions):
    low, high = 0, len(regions)
    while low < high:
        mid = (low + high) // 2
        region = regions[mid]
        if pos[0][-1] < region[0]:
            # position is located on a chromosome/o-genome 
            # with lesser id than current region
            high = mid
        elif pos[0][-1] > region[0]:
            # position is located on a chromosome/o-genome 
            # with higher id than current region
            low = mid + 1
        else:
            # position is located on same genome as region ...
            if pos[1] < region[2]:
                # ... but upstream of region
                high = mid
            elif pos[1] > region[3]:
                # ... but downstream of region
                low = mid + 1
            else:
                # and is also located within region
                return region
    return None

###
def remove_intragenic_snps(snp_d, intragenic_regions):
    for pos in snp_d:
        region = find_region(pos, intragenic_regions)
        if region is None:
            region = ['None']
        print '%c %i' % pos, ' '.join(map(str, region))
    
    return snp_d

###
def main(argv):

    samfile = pysam.Samfile(argv[0], 'rb')
    fo = sys.stdout
    fo = open('%s.covered_genes.csv' % argv[0].rstrip('.bam'), 'w')    
    # fo2 = open('%s.covered_genes_with_reads.csv' % argv[0].rstrip('.bam'), 'w')
    fo2=sys.stdout

    # fo.write('%s\n' % ','.join(COL_HEADERS))
    
    # tair10_genes.gff
    intragenic_regions = gff_helpers.read_intragenic_regions(open(argv[1]))
    # ped-0-snps_no-indels.txt
    snp_d = SNPDict(open(argv[2]))  
    
    snp_d = remove_intragenic_snps(snp_d, intragenic_regions)
    
    
    # tair10_genes.gff
    # process_gff(open(argv[3]), polymorphs, snp_d, samfile, fo, fo2, min_reads=MIN_NREADS)
    # fo2.close()
    
    
    fo.close()
    samfile.close()
    return None



if __name__ == '__main__': main(sys.argv[1:])
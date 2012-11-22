#!/usr/bin/env python
'''
Created on Nov 16, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


def main(argv):
    
    genome_path = '/home/schudoma/projects/ngs/tair10/TAIR10_chr%c.fas'
    snpfile = '/home/schudoma/projects/ngs/Ath_Ped-0-Col-0-graft_RNA-Seq_Scheible/ath_ped-0_snps_no-N_no-indel.txt'
    
    current_contig = None
    refseq = None
    
    for line in open(snpfile):
        line = line.strip().split('\t')
        # Ped-0   1       110     G       T       27      3       1       1
        contig, pos, refbase, mutation = line[1], int(line[2]), line[3], line[4]
        if contig == '6': 
            contig = 'C'
        elif contig == '7':
            contig = 'M'
        if refseq is None or current_contig != contig:
            refseq = ''.join(open(genome_path % contig).readlines()[1:]).replace('\n', '')
            current_contig = contig
        # if refseq[pos] == refbase, then the offset used in the snpfile is 0
        if refseq[pos] != refbase:
            print 'First error at %s:%i (%c, should be: %c)' % (contig, pos, refseq[pos], refbase),
            # it's not 0, testing offset=1
            if refseq[pos - 1] != refbase:
                print
                print 'Second error at %s:%i (%c, should be: %c)' % (contig, pos - 1, refseq[pos - 1], refbase)
                break
            else:
                print '... but offset=1 seems fine: %c = %c' %  (refseq[pos - 1], refbase)
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

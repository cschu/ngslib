#!/usr/bin/env python
'''
Created on Mar 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


def main(argv):
    N_ACCESSIONS = 80.0
    CUTOFF = 0.05
    OUT = sys.stdout
    
    snpcount = 0
    for line in open(argv[0]):
        # if snpcount >= 10: break
        line = line.strip().split()
        (contig, pos, col), accessions = line[:3], line[2:]
        ped = accessions[66]
        
        if col != ped:
            # we need col and ped to differ
            if col in 'ACGT' and ped in 'ACGT':
                # we need col and ped to carry substitution SNPs only
                col_count = accessions.count(col)
                ped_count = accessions.count(ped)
                # at least 5% of accessions have to support either SNP 
                if (col_count / N_ACCESSIONS) >= CUTOFF and (ped_count / N_ACCESSIONS) >= CUTOFF:
                    line += map(str, [col_count, ped_count, col, ped])
                    OUT.write(' '.join(line) + '\n')
                    snpcount += 1
    OUT.write('%i multicovered SNPs found.\n' % snpcount)
        
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

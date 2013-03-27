#!/usr/bin/env python
'''
Created on Mar 19, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

from collections import Counter

def main(argv):
    
    n_substitution_snps = 0
    n_col_singletons = 0
    n_singletons = 0
    
    for line in open(argv[0]):
        line = line.strip().split('\t')
        
        contig, pos, col_snp = line[:3]
        
        
        counts = Counter(line[2:])
        if len(set(counts.elements()) - set(['A', 'C', 'G', 'T', 'U'])) > 0:
            continue
        n_substitution_snps += 1
        other_singletons = [x for x in counts.items() if x[1] == 1 and x[0] != col_snp]
        
        if counts[col_snp] == 1:
            n_col_singletons += 1
            
        if len(other_singletons) > 0:
            n_singletons += len(other_singletons)
        
        
        
        
        
    print n_col_singletons, 'Col-0 singletons and',
    print n_singletons, 'other singletons',      
    print 'in', n_substitution_snps, '(%.3f)' % (float(n_col_singletons)/n_substitution_snps)
            
        
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

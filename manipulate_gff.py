#!/usr/bin/env python
'''
Created on Nov 15, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re

# Chr1    .       snp     4879    4880    .       .       .       quality=38;supporting_reads=25;concordance=0.925926;avg_align_overl_reads=1.0;refbase=C;mutation=T
# Chr1    TAIR10  gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010    2

def main(argv):
    
    fn = argv[0]
    agi = re.compile('ID=AT[12345CM]G[0-9]{5}')
    out = sys.stdout
    
    for line in open(fn):
        # print line
        line = line.strip().split('\t')
        newline = '\t'.join(line[:8])
        attributes = line[8]
        # print line[17]
        new_attributes = ['gene_start=%s' % line[12],
                          'gene_end=%s' % line[13],
                          'gene_strand=%s' % line[15],
                          'gene_' + agi.search(line[17]).group()]
        attributes += ';' + ';'.join(new_attributes)
        newline += '\t' + attributes
        out.write(newline + '\n')
        
    pass

if __name__ == '__main__': main(sys.argv[1:])

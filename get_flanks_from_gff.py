#!/usr/bin/env python
'''
Created on Jan 24, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

from get_introns_from_gff import get_attributes

def main(argv):
    
    out = sys.stdout
    
    cds_d = {}
    try:
	flanksize = int(argv[1]) 
    except:
	flanksize = 300

    
    for line in open(argv[0]):
        line = line.strip().split('\t')
        attributes = get_attributes(line[8])
        
        if line[2] == 'CDS':
            id_ = attributes['Parent'].split(',')[0]
            cds_d[id_] = cds_d.get(id_, []) + [int(line[3]), int(line[4])]
    
    for id_, coding_sequences in sorted(cds_d.items()):
        p5, p3 = min(coding_sequences), max(coding_sequences)
        
        line = ['Chr%c' % id_[2], 'TAIR10', 'five_prime_flank', str(max(1, p5 - flanksize)), 
                str(p5), '.', '.', '.', 'Parent=%s' % id_]
        out.write('%s\n' % '\t'.join(line))
        line = ['Chr%c' % id_[2], 'TAIR10', 'three_prime_flank', str(p3), 
                str(p3 + flanksize), '.', '.', '.', 'Parent=%s' % id_]
        out.write('%s\n' % '\t'.join(line))
        
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

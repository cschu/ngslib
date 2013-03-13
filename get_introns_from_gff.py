#!/usr/bin/env python
'''
Created on Jan 24, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

def get_attributes(attributes):
    return dict([tuple(field.split('=')) for field in attributes.rstrip(';').split(';')])


def main(argv):
    
    out = sys.stdout
    
    mrna = None
    last_exon = None
    for line in open(argv[0]):
        line = line.strip().split('\t')
        
        
        attributes = get_attributes(line[8])
        
        if line[2] != 'exon':            
            # print line, 'not exon'
            mrna = tuple(map(int, line[3:5])) + (attributes['ID'],)
            last_exon = None
        else:
            current_exon = tuple(map(int, line[3:5])) + (attributes['Parent'],)
            if not last_exon is None: 
                line[2] = 'intron'
                if line[6] == '+':
                    line[3] = str(last_exon[1] + 1)
                    line[4] = str(current_exon[0] - 1)
                else:
                    line[4] = str(last_exon[0] - 1)
                    line[3] = str(current_exon[1] + 1)
                if mrna[2] == last_exon[2] and mrna[2] == current_exon[2]:
                    out.write('%s\n' % '\t'.join(line))
                else:
                    out.write('###Cannot write intron: mrna=%s last=%s current=%s\n' % \
                              tuple(map(str, [mrna, last_exon, current_exon])))
                    # sys.exit(1)
            last_exon = current_exon
                
            
              
             
        
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

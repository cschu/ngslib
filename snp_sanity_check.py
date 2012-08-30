#!/usr/bin/env python
'''
Created on Aug 21, 2012

@author: Chris
'''

import sys

def get_gffstring(contig, pos):
    return '\t'.join([contig, '.', 'snp', pos, pos, '.', '.', '.', '.'])
    

def main(argv):
    
    header = None
    for line in open(argv[0]):
        if header is None:
            header = line
        elif line.startswith('#'):
            continue
        else:
            line = line.split(';')
            print get_gffstring('Chr%c' % line[0], line[1])
    
    
    
    pass

    
if __name__ == '__main__': main(sys.argv[1:])
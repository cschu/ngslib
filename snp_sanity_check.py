#!/usr/bin/env python
'''
Created on Aug 21, 2012

@author: Chris
'''

import sys

def get_gffstring(contig, pos):
    return '\t'.join([contig, '.', 'snp', pos, pos, '.', '.', '.', '.'])
    

def main(argv):
    
    for line in open(argv[0]):
        line = line.split(';')
        print get_gffstring('Chr%c' % line[0], line[1])
    
    
    
    pass

    
if __name__ == '__main__': main(sys.argv[1:])
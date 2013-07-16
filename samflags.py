#!/usr/bin/env python
'''
Created on Jun 4, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


SAMFLAGS = [(1, 'multiple segments'),
            (2, 'each segment properly aligned'),
            (4, 'segment unmapped'),
            (8, 'next segment unmapped'),
            (16, 'SEQ is reverse complemented'),
            (32, 'SEQ of next segment is reverse complemented'),
            (64, 'first segment'),
            (128, 'last segment'),
            (256, 'secondary alignment'),
            (512, 'qc failed'),
            (1024, 'PCR/optical duplicate')] 


def main(argv):
    
    flags = int(argv[0])
    
    
    for flag, desc in SAMFLAGS:
        if (flags & flag) == flag:
            print flag, desc
        
        
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

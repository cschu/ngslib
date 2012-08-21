#!/usr/bin/env python

'''
Created on Aug 21, 2012

@author: Chris
'''

import sys

import pysam


def main(argv):
    
    samfile = pysam.Samfile(argv[0], 'rb')
    
    
    for col in samfile.pileup('Chr1', start=4880, end=4881, stepper='all'):
        
        if col.pos == 4880:                                              
            for read in col.pileups:
                ali = read.alignment
                print dict(ali.tags)['RG'], ali.pos, ali.alen, ali.rlen, ali.qlen
                
    
    samfile.close()
            
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
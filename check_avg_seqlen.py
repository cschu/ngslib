#!/usr/bin/python
'''
Created on Jul 11, 2012

@author: schudoma
'''
import sys

import pysam


def main(argv):
    samfile = pysam.Samfile(argv[0], 'rb')
    
    lengths = [read.rlen for read in samfile.fetch()]
    print argv[0], sum(lengths)/float(len(lengths))
    
    
    samfile.close()
    
    pass


if __name__ == '__main__': main(sys.argv[1:])
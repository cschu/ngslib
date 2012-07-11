#!/usr/bin/python
'''
Created on Jun 28, 2012

@author: schudoma
'''

import sys
import pysam


def main(argv):
    reads = pysam.Samfile(argv[0])
    
    for read in reads:
        if read.is_proper_pair and read.is_read1:
            pos = reads.tell()
            
            try:
                mate = reads.mate(read)
            except ValueError:
                continue
            finally:
                reads.seek(pos)
            
            print read.qname, read.tid, 'read1'
            print mate.qname, mate.tid, 'read2'
            print
            print read
            print mate
            sys.exit()
    
    

if __name__ == '__main__': main(sys.argv[1:])
          

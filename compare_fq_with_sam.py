#!/usr/bin/env python
'''
Created on Jun 18, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

def main(argv):
    
    fqids = set()
    for line in open(argv[0]):
        if line.startswith('@HWI'):
            fqids.add(line.strip().strip('@').split()[0])
    
    
    
    for samfile in argv[1:]:        
        sam_reads = set([])
        for line in open(samfile):
            line = line.strip().split()
            sam_reads.add(line[0])
        
        missing_reads = len(fqids.difference(sam_reads))        
        print samfile, missing_reads, len(fqids), str(missing_reads == 0).upper() 
        
    pass

if __name__ == '__main__': main(sys.argv[1:])

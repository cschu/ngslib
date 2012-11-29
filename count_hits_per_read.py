#!/usr/bin/env python
'''
Created on Nov 22, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import pysam

def main(argv):
    bamfile = pysam.Samfile(argv[0], 'rb')
    counts = {}
    for alignedRead in bamfile.fetch():
        if alignedRead.flag & 0x40 == 0x40 or alignedRead.flag == 0:
            counts[alignedRead.qname] = counts.get(alignedRead.qname, 0) + 1
    bamfile.close()
    
    out = open(argv[0] + '.mult_counts', 'wb')
    for item in sorted(counts.items()):
        out.write('%s\t%i\n' % item)
    out.close()
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

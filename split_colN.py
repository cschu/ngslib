#!/usr/bin/env python
'''
Created on Nov 29, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import pysam

def main(argv):
    bamfile = pysam.Samfile(argv[0], 'rb')    
    reads = {}
    for alignedRead in bamfile.fetch():
        reads[alignedRead.qname] = reads.get(alignedRead.qname, []) + [alignedRead]
    bamfile.close()
    
    total = sum([len(v) for v in reads.values()])
    splits = 5
    splitsize = float(total)/splits
    

    bamout = None
    cursize = 0
    i = 0
    for k, v in sorted(reads.items()):
        if bamout is None or cursize > splitsize:
            if bamout is not None:
                bamout.close
            bamout = pysam.Samfile(argv[0] + '.%02i' % i, 'wb')
            cursize = 0
            i += 1
        for alignedRead in v:
            bamout.write(alignedRead)
            cursize += 1
    try:
        bamout.close()
    except:
        pass
    
    
    
        
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

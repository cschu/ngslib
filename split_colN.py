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
    sys.stdout.write('Gathering read data...\n')
    sys.stdout.flush()
    for alignedRead in bamfile.fetch(until_eof=True):
        reads[alignedRead.qname] = reads.get(alignedRead.qname, []) + [alignedRead]

    header = bamfile.header
    bamfile.close()
    
    total = sum([len(v) for v in reads.values()])
    splits = 5
    splitsize = float(total)/splits
    
    bamout = None
    cursize = 0
    i = 0
    written = 0
    for k, v in sorted(reads.items()):
        if bamout is None or cursize > splitsize:
            if bamout is not None:
                bamout.close
            
            sys.stdout.write('Reads in last file: %i\n' % cursize)
            sys.stdout.write('Total reads written: %i, left: %i\n' % \
                                 (written, total-written))
            sys.stdout.write('Now writing reads to file %s.%02i\n' % \
                                 (argv[0], i))
            sys.stdout.flush()
            bamout = pysam.Samfile(argv[0] + '.%02i' % i, 'wb', header=header)
            cursize = 0
            i += 1
        for alignedRead in v:
            bamout.write(alignedRead)
            cursize += 1
            written += 1
    
    sys.stdout.write('Reads in last file: %i\n' % cursize)
    sys.stdout.write('Total reads written: %i, left: %i\n' % \
                         (written, total-written))
            
    try:
        bamout.close()
    except:
        pass
    
    
    
        
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

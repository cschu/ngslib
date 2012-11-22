#!/usr/bin/env python
'''
Created on Nov 15, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


def read_gff(lines):
    return set([tuple(line.strip().split('\t')) for line in lines.split('\n') if len(line) > 0])

def main(argv):
    gffset = read_gff(open(argv[0]).read())
    for fn in argv[1:]:        
        gffset = gffset.union(read_gff(open(fn).read()))
    for item in gffset:
        if len(item) != 9:
            print item
            sys.exit()
    
    for line in sorted(list(gffset), key=lambda x:(x[0], int(x[3]), int(x[4]))):
        sys.stdout.write('\t'.join(line) + '\n')
    
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

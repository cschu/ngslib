#!/usr/bin/env python
'''
Created on Aug 30, 2012

@author: Chris
'''
import sys

def fq2fasta(fq_in, fas_out=sys.stdout):
    while True:
        try:
            header = fq_in.readline().strip().split()[0]
            seq = fq_in.readline().strip()
            dummy = fq_in.readline().strip()
            qual = fq_in.readline().strip()
        except:
            break
        out_str = '>%s|%s\n%s\n' % (header, qual, seq)
        out_str = '>%s\n%s\n' % (header, seq)
        fas_out.write(out_str)
        pass
    pass


def main(argv):
    for arg in argv:
        fq2fasta(open(arg, 'r'))
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python
'''
Created on Apr 5, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


# 1       1       C       Z



def main(argv):
    LINEWIDTH = 79
    
    contig = None
    colfile = open('col_genome.fa', 'w')
    pedfile = open('ped_genome.fa', 'w')
    
    colbuffer = ''
    pedbuffer = ''
    for line in open(argv[0]):
        line = line.strip().split('\t')
        if contig is None or line[0] != contig:
            contig = line[0]
            header = '>Chr%c' % contig
            colfile.write(header + '\n')
            pedfile.write(header + '\n')
        if len(colbuffer) == LINEWIDTH:
            colfile.write(colbuffer + '\n')
            pedfile.write(pedbuffer + '\n')
            colbuffer, pedbuffer = '', ''
        colbuffer += line[2]
        pedbuffer += line[3]
    if len(colbuffer) > 0:
        colfile.write(colbuffer + '\n')
        pedfile.write(pedbuffer + '\n')
    
    colfile.close()
    pedfile.close()
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

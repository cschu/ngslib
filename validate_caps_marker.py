#!/usr/bin/env python
'''
Created on Mar 19, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


def main(argv):
    
    sequences = {}
    for line in open(argv[0]):
        line = line.strip()
        print line, line.startswith('>'), line[0]
        if line.startswith('>'):
            key = line[1:]
        else:
            sequences[key] = line
    for k, v in sequences.items(): 
        print '>' + k
        start = int(k.split('_')[1])
        while True:
            p = v.find('(') 
            if p == -1: break
            print 'SNP at', start+p, v[p:p+5]            
            #v = v[:p] + v[p+1] + v[p+5:]
            v = v[:p] + 'X' + v[p+5:]
            print v, len(v)
            # break
            
                
        
        # print v
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

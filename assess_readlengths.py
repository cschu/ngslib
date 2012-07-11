#!/usr/bin/python

import os
import sys
import math

###
def do_something():
    return None 

###
def main(argv):
    sequences = []
    state = None
    for line in open(argv[0]):
        # print state
        if line.startswith('@'):
            state = 'header'
        elif state == 'header':
            sequences.append(line.strip())
            print sequences[-1]
            state = 'sequence'
        elif state == 'sequence':
            state = '+'
        elif state == '+':
            state = 'quality'        
        pass
    
    print sequences

    return None

if __name__ == '__main__': main(sys.argv[1:])

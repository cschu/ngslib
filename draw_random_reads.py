#!/usr/bin/env python
'''
Created on Aug 22, 2012

@author: schudoma
'''

import sys
import random

def main(argv):
    
    data = open(argv[0]).readlines()
    rnd = range(len(data))
    random.shuffle(rnd)
    N = int(argv[1])
    #rnd = rnd[:N]
    for i in rnd[:N]:
        print data[i]
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
#!/usr/bin/env python
'''
Created on Apr 17, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

import random
from collections import Counter


def main(argv):
    
    N = 30000
    N_SAMPLES = 1000
    urn = range(N)
    s = 300
    
    counts = []
    for i in xrange(N_SAMPLES):
        sample = [random.choice(urn) for j in xrange(s)]
        c = Counter(sample)
        c.subtract(urn)
        counts.append(len(list(c.elements())))
    
    counts = map(lambda x: x/float(s), counts)
    print sum(counts) / N_SAMPLES
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

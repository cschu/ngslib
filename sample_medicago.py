#!/usr/bin/env python
'''
Created on Jun 26, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import random


def main(argv):
    
    fn1 = argv[0]
    fn2 = argv[1]
    
    print 'Counting lines...'
    n_items = sum([1 for line in open(fn1)]) / 4
    print 'Making item vector...'
    rnd_vector = range(n_items)
    print 'Shuffling item vector...'
    random.shuffle(rnd_vector)
    print 'Sorting reduced item vector...'
    rnd_vector = sorted(rnd_vector[:int(argv[2])])
    
    l_count = 0
    fi1 = open(fn1)
    fi2 = open(fn2)
    
    fo1 = open(fn1.replace('.fastq', '.sample_%04i.fastq' % int(argv[3])), 'wb')
    fo2 = open(fn2.replace('.fastq', '.sample_%04i.fastq' % int(argv[3])), 'wb')
    
    print 'Scanning files...'
    while len(rnd_vector) > 0:
        lines1 = [fi1.readline(), fi1.readline(), fi1.readline(), fi1.readline()]
        lines2 = [fi2.readline(), fi2.readline(), fi2.readline(), fi2.readline()]
        if l_count == rnd_vector[0]:
            print l_count, len(rnd_vector), lines1[0], lines2[0]
            fo1.write(''.join(lines1))
            fo2.write(''.join(lines2))
            fo1.flush()
            fo2.flush()
            del rnd_vector[0]
        l_count += 1
    fo1.close()
    fo2.close()
    
    fi1.close()
    fi2.close()
    
        
        
    
    
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

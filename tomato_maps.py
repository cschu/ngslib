#!/usr/bin/env python
'''
Created on Sep 3, 2012

@author: schudoma
'''
import os
import sys
import glob

def main(argv):
    
    mapping = {}
    
    files = glob.glob('*.map')
    for fn in files:
        stat = os.stat(fn)
        nlines = len(open(fn).readlines())
        fn = fn.split('_')
        set1, set2 = '_'.join(fn[:3]), '_'.join(fn[3:]).strip('.map')
        mapping[(set1, set2)] = (stat.st_size, nlines)
        
        
        
    for k, v in sorted(mapping.items()):
        print '%s,%s,%i,%i' % (k[0], k[1], v[0], v[1])
    
    
    # CTS2_snRNA_R2.unpaired_TSS1_snRNA_R1_up.rev.CTS.map
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
    
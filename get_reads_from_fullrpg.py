#!/usr/bin/env python
'''
Created on Mar 1, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import pickle

def main(argv):
    
    # agi = argv[0]
    rpg_files = argv[0:]
    
    for agi in ['AT4G20270', 'AT3G56040']:
        for fn in rpg_files:
            dat = pickle.load(open(fn))
        
            ped = [x[0] for x in dat[agi] if x[2]]
            col = [x[0] for x in dat[agi] if not x[2]]
        
            conflict = set(ped).intersection(set(col))
            ped_uniq = sorted(set(ped) - conflict)
            col_uniq = sorted(set(col) - conflict)
            
            is_paired = str(int(not 'unpaired' in fn))
        
            for rid in ped_uniq:
                sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % (agi, rid, 'ped', str(ped.count(rid)), is_paired))
            for rid in col_uniq:
                sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % (agi, rid, 'col', str(col.count(rid)), is_paired))
    
    #AT4G20270       HWI-ST377:127:D0PHGACXX:7:1105:21207:11467      ped     2       0

    
    pass

if __name__ == '__main__': main(sys.argv[1:])

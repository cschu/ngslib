#!/usr/bin/env python
'''
Created on Dec 21, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)

This script searches cucumber-arabidopsis orthologs. 
'''
import sys
import math
# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
# Csa1M074400.1   AT5G46690.1     93.33   30      2       0       51      80      276     305     0.001   44.1

def read_tab_blast(open_fn):
    
    csa_mobiles = {}
    for line in open_fn:
        if not line.startswith('#'):
            line = line.strip().split('\t')
            csa, ath = line[:2]            
            csa_mobiles[csa] = csa_mobiles.get(csa, []) + ([[ath] + map(float, line[2:])])
    return csa_mobiles


def find_orthologues(csa, ath_data, n=3):
    # print ath_data
    maxlen = max([x[2] for x in ath_data])
    # print maxlen
    candidates = []
    for i, candidate in enumerate(ath_data):
        if candidate[2] == maxlen:
            candidates.append(candidate)
            continue
        # e-value <= 10e-3
        if candidate[9] > 0.001:
            # print 'evalue', candidate 
            continue
        # identity >= 80%
        if round(candidate[1]) < 80.0:
            # print 'identity', candidate
            continue
        # alignment length >= 0.8 * max_alignment_length
        if candidate[2] < 0.8 * maxlen:
            # print 'length', candidate
            continue
        candidates.append(candidate)
    return candidates

    
def main(argv):
    
    data = read_tab_blast(open(argv[0]))
    for k, v in data.items():
        candidates = find_orthologues(k, v)
        # print k
        # print v
        print '#', k, len(v), len(candidates), float(len(candidates))/len(v)
        # print candidates
        for c in candidates:
            print k, '\t'.join(map(str, c)) 
        # break
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

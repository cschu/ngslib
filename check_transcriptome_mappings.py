#!/usr/bin/env python
'''
Created on Jun 18, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re
import os

from collections import defaultdict


def main(argv):
    
    mapping = set([])
    not_mapping = set([])
    transcripts = defaultdict(set)
    
    try:
        sample_id = os.path.basename(argv[0])[:4]
    except:
        sample_id = 'DUMMY'
    
    for fn in argv:        
        for line in open(fn):
            if line.startswith('@'):
                continue
            fields = line.strip().split()
            # print fields
            is_mapped = (int(fields[1]) & (4 == 0)) and re.search('NM:i:[01]', line)        
            mate_id = fields[0].strip()
        
            if is_mapped:
                mapping.add(mate_id)
                tid = fields[2].strip()
                transcripts[tid].add(mate_id)
            else:
                not_mapping.add(mate_id)
                
    pair_mismatches = mapping.intersection(not_mapping)
    mapping = mapping.difference(pair_mismatches)
    
    for tid in transcripts:
        transcripts[tid] = transcripts[tid].difference(pair_mismatches)
    
    out = open('%s_HOST_MAPPING.txt' % sample_id, 'wb')
    for rid in sorted(list(mapping)):
        out.write('%s\n' % rid)
    out.close()
    out = open('%s_HOST_NOTMAPPING.txt' % sample_id, 'wb')
    for rid in sorted(not_mapping):
        out.write('%s\n' % rid)
    out.close()
    out = open('%s_TRANSCRIPTS.txt' % sample_id, 'wb')
    for tid in transcripts:
        n_reads = len(transcripts[tid])
        if n_reads > 0:
            out.write('%s\t%i\n' % (tid, n_reads))            
    out.close()
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python
'''
Created on Jun 20, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os
import time, datetime

from collections import defaultdict

from trie import DuplicateReadFinder

def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')

""" CCL1_unpaired_all.fastq
CCL1_R1_paired.fastq    CCL1_R2_paired.fastq
"""

def read_seqs(fn, drf, wanted):
    key = None    
    for line in open(fn):
        if line.startswith('@HWI'):
            key = line.strip().strip('@').split()[0]
        elif key is not None:
            if key in wanted:
                drf.add_sequence(wanted[key], line.strip())            
            key = None
    pass

def main(argv):
    
    fn_list = argv[0]
    seqdir = argv[1]
    seqfile = open(argv[2] + '.fa', 'wb')
    
    drf = DuplicateReadFinder()
    wanted = dict()
    for line in open(fn_list):
        line = line.strip()
        id_, sources = line.split('\t')
        sources = sources.split(',')
        if len(sources) == 2:
            kk = 'p'
        elif sources[0].startswith('R'):
            kk = 'u'
        else:
            kk = ''
        for k in sources:
            wanted[id_] = '%s_%s%s' % (id_, k, kk)
    
    
    setname = seqdir.split('_')[1]    
    sys.stderr.write('%s> Reading sequences and checking for duplicates...\n' % get_timestamp())
    
    for fn in ['%s_unpaired_all.fastq', '%s_R1_paired.fastq', '%s_R2_paired.fastq']:
        read_seqs(os.path.join(seqdir, fn % setname), drf, wanted)
            
    sys.stderr.write('%s> Writing sequences...\n' % get_timestamp())
    drf.write_sequences(seqfile, modify_output=lambda x:x.replace('U', 'T'))                
    seqfile.close()
    sys.stderr.write('%s> Done.\n' % get_timestamp())
        
        


    pass

if __name__ == '__main__': main(sys.argv[1:])

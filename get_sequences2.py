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

def read_seqs(fn):
    key = None
    seqs = {}
    for line in open(fn):
        if line.startswith('@HWI'):
            key = line.strip().strip('@').split()[0]
        elif key is not None:
            seqs[key] = line.strip()
            key = None
    return seqs


def main(argv):
    
    fn_list = argv[0]
    seqdir = argv[1]
    seqfile = open(argv[2] + '.fa', 'wb')

    setname = seqdir.split('_')[1]
    sys.stderr.write('%s> Reading sets...\n' % get_timestamp())
    seq_d = {'U': read_seqs(os.path.join(seqdir, 
                                         '%s_unpaired_all.fastq' % setname)),
             'R1': read_seqs(os.path.join(seqdir,
                                          '%s_R1_paired.fastq' % setname)),
             'R2': read_seqs(os.path.join(seqdir,
                                          '%s_R2_paired.fastq' % setname))}
    sys.stderr.write('%s> Checking for duplicates...\n' % get_timestamp())
    # seqfile = open(os.path.basename(fn_list) + '.fa', 'wb')        
    
    drf = DuplicateReadFinder()    
    for line in open(fn_list):
        line = line.strip()
        id_, sources = line.split('\t')
        sources = sources.split(',')
        for k in sources:
            if len(sources) == 2:
                kk = 'p'
            elif sources[0].startswith('R'):
                kk = 'u'
            else:
                kk = ''
            if id_ in seq_d[k]:
                drf.add_sequence(id_ + '_' + k + kk, seq_d[k][id_])
                # seqfile.write('>%s_%s\n%s\n' % (id_, k, seq_d[k][id_]))
    
    sys.stderr.write('%s> Writing sequences...\n' % get_timestamp())
    drf.write_sequences(seqfile)                
    seqfile.close()
    sys.stderr.write('%s> Done.\n' % get_timestamp())
        
        


    pass

if __name__ == '__main__': main(sys.argv[1:])

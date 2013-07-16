#!/usr/bin/env python
'''
Created on Jun 20, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os


from collections import defaultdict

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

    setname = seqdir.split('_')[1]
    print 'Reading sets...'
    seq_d = {'U': read_seqs(os.path.join(seqdir, 
                                         '%s_unpaired_all.fastq' % setname)),
             'R1': read_seqs(os.path.join(seqdir,
                                          '%s_R1_paired.fastq' % setname)),
             'R2': read_seqs(os.path.join(seqdir,
                                          '%s_R2_paired.fastq' % setname))}
    print 'Generating seqfile...'
    seqfile = open(os.path.basename(fn_list) + '.fastq', 'wb')
    for line in open(fn_list):
        line = line.strip().split()
        for k in ['R1', 'R2', 'U']:
            if line[0] in seq_d[k]:
                seqfile.write('>%s_%s\n%s\n' % (line[0], k, seq_d[k][line[0]]))
    seqfile.close()
        
        


    pass

if __name__ == '__main__': main(sys.argv[1:])

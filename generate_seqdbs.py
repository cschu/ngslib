#!/usr/bin/env python
'''
Created on Jul 15, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import os
import sys
import anydbm
import time, datetime

def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
    
def read_seqs(fn, dbm):
    key = None    
    for line in open(fn):
        if line.startswith('@HWI'):
            key = line.strip().strip('@').split()[0]
        elif key is not None:            
            dbm[key] = line.strip()
            key = None
    dbm.close()
    pass

def main(argv):
    
    seqdir = argv[0]
    setname = seqdir.split('_')[1]
    
    sys.stderr.write('%s> Reading sets...\n' % get_timestamp())
    read_seqs(os.path.join(seqdir, '%s_unpaired_all.fastq' % setname), 
              anydbm.open(argv[1] + '.U.dat', 'c'))    
    read_seqs(os.path.join(seqdir, '%s_R1_paired.fastq' % setname),
              anydbm.open(argv[1] + '.R1.dat', 'c'))    
    read_seqs(os.path.join(seqdir, '%s_R2_paired.fastq' % setname),
              anydbm.open(argv[1] + '.R2.dat', 'c'))
    sys.stderr.write('%s> Done.\n' % get_timestamp())
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

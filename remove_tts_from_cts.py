#!/usr/bin/env python

'''
Created on Sep 5, 2012

@author: schudoma
'''
import sys
import os
import glob
import subprocess as subp


def make_seqdict(path, pattern=''):
    print 'Generating seqdict...'
    seqdict = {}
    for fn in glob.glob(os.path.join(path, pattern)):
        # if '_R2.paired' in fn:
        #    continue        
        print 'Processing', fn
        for line in open(fn):
            if line.startswith('>'):
                line = line.strip()[1:]
                #if line in seqdict:                    
                #    print 'Problem:', line, 'is not unique.', seqdict[line]
                #    sys.exit(1)                                            
                seqdict[line] = seqdict.get(line, []) + [os.path.basename(fn)]
                # seqdict[line] = os.path.basename(fn)                
    return seqdict

def get_read_ids(path, col):
    print 'Getting read ids from', path, '...'
    read_ids = set([])
    for fn in glob.glob(os.path.join(path, '*.map')):
        print 'Processing', fn
        for line in open(os.path.join(path, fn)):
            line = line.strip().split('\t')
            read_ids.add(line[col])
    return list(read_ids)

def read_fasta(open_fn):
    seqs = {}
    last = None
    for line in open_fn:
        line = line.strip()
        if line.startswith('>'):
            last = line[1:]
        else:
            seqs[last] = line
    return seqs

def write_sequences(seqdict, seqpath, out=sys.stdout):
    
    file_d = {}
    print 'Generating filename -> read id map... (%i items)' % len(seqdict)
    # items = sorted([(v, k) for k, v in seqdict.items()])
    items = []
    for k, v in seqdict.items():
        items.extend([(v_i, k) for v_i in v])
    
    current = None
    buffer_ = []
    
    for item in sorted(items):
        # print item
        if current is None or item[0] != current:
            print 'New file'
            if not current is None:
                print 'Dumping buffer'
                file_d[current] = buffer_
            current = item[0]
            buffer_ = []
        buffer_.append(item[1])
    print 'Dumping buffer'        
    file_d[current] = buffer_
        
    # print len(seqdict)
    print sum(map(len, seqdict.values()))
    print sum(map(len, file_d.values()))
    
    # for k, v in seqdict.items():
    #    file_d[v] = file_d.get(v, []) + [k]
    # print 
    for fn, rids in sorted(file_d.items()):
        seqs = read_fasta(open(os.path.join(seqpath, fn)))
        for rid in rids:
            out.write('>%s|%s\n%s\n' % (rid, fn.rstrip('.fas'), seqs[rid]))
    pass


def main(argv):
    
    ngspath = '/home/schudoma/ngstemp/cschudoma'
    workpath = os.path.join(ngspath, 'Sly-Cus_parasitic_RNA-Seq_Kragler')    
    tts_cts = os.path.join(workpath, 'bowtie-TTS-CTS-results')
    cts_tts = os.path.join(workpath, 'bowtie-CTS-TTS-results')
    readpath = os.path.join(workpath, 'fasta')
    
    # step 0: get all sequences from CTS reads
    seqdict = make_seqdict(readpath, pattern='CTS*.fas')
    # print max(map(len, seqdict.values())) # -- debug
    
    # step 1: get all reads that mapped in CTS->TTS or TTS->CTS
    unique_read_ids = list(set(get_read_ids(tts_cts, 2) + get_read_ids(cts_tts, 0)))
    
    # step 2: remove all reads obtained in step 1 from the set of all CTS reads (from step 0)
    for rid in unique_read_ids:
        if rid in seqdict:
            del seqdict[rid]
        else:
            print 'Weird:', rid, 'is not contained in seqdict.'
            
    # step 3: write all remaining reads (this is the set CTS\(CTS^TTS))            
    write_sequences(seqdict, readpath, out=open('CTS_not_in_TTS.fas', 'wb'))
    
            
    
    
     
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
    
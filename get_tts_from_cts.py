#!/usr/bin/env python
'''
Created on Sep 6, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os



from remove_tts_from_cts import make_seqdict, get_read_ids, read_fasta #, write_sequences

def read_fq(open_fn):
    seqs = {}
    last = None
    for line in open_fn:
        line = line.strip()
        if line.startswith('@HWI'):
            
            # the .split()[0] results from Casava 1.8 format
            last = line.strip().split()[0]
            seqs[last] = ''
        else:
            seqs[last] += line.strip() + '\n'
    return seqs


def write_sequences(seqdict, seqpath, out=sys.stdout):
    
    print seqdict
    
    file_d = {}
    print 'Generating filename -> read id map... (%i items)' % len(seqdict)    
    items = []
    for k, v in seqdict.items():
        items.extend([(v_i, k) for v_i in v])
    
    current = None
    buffer_ = []
    
    for item in sorted(items):
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
        
    print sum(map(len, seqdict.values()))
    print sum(map(len, file_d.values()))
    
    # print file_d.items()
    # sys.exit(1)
    
    for fn, rids in sorted(file_d.items()):   
        fqfn = fn.replace('.fas', '.fq')     
        print seqpath, fqfn
        
        seqs = read_fq(open(os.path.join(seqpath, fqfn)))
        # print seqs.keys()[:10]
        # sys.exit(1)
        out = open('%s_in_TSS.fq' % fn.rstrip('.fas'), 'wb')
        out2 = None
        seqs2 = {}                
        if '.paired.fq' in fqfn:
            fn2 = fqfn.replace('R1', 'R2') #.rstrip('.fas')
            out2 = open('%s_in_TSS.fq' % fn2.rstrip('.fas'), 'wb')
            seqs2 = read_fq(open(os.path.join(seqpath, fn2.replace('.fas', '.fq'))))
        
        for rid in rids:
            out.write('%s\n%s\n' % (rid, seqs[rid].strip()))
            if out2:
                out2.write('%s\n%s\n' % (rid, seqs2[rid].strip()))
        pass
            
    pass

def main(argv):
    
    ngspath = '/home/schudoma/ngstemp/cschudoma'
    workpath = os.path.join(ngspath, 'Sly-Cus_parasitic_RNA-Seq_Kragler')
    cts_tts = os.path.join(workpath, 'bowtie-CTS-TTS-results')
    readpath = os.path.join(workpath, 'fasta')
    fqpath = os.path.join(workpath, 'fastq')
    
    # step 0: get all sequences from CTS reads
    seqdict = make_seqdict(readpath, pattern='CTS*.fas')
    
    # step 1: get all reads that mapped in CTS->TTS 
    unique_read_ids = list(set(get_read_ids(cts_tts, 0)))
    
    # step 2: remove all reads obtained in step 1 from the set of all CTS reads (from step 0)
    keep_seq = {}
    
    print 'Gathering sequences for export...'
    for rid in unique_read_ids:
        if rid in seqdict:
            keep_seq[rid] = seqdict[rid]
            # del seqdict[rid]
        else:
            print 'Weird:', rid, 'is not contained in seqdict.'
    
    # del seqdict
    print 'keep_seq has %i entries' % len(keep_seq)       
    # step 3: write all remaining reads (this is the set CTS\(CTS^TTS))            
    write_sequences(keep_seq, fqpath)

if __name__ == '__main__': main(sys.argv[1:])

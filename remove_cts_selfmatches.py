#!/usr/bin/env python
'''
Created on Sep 6, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import os
import sys
import glob

from collections import defaultdict

def make_seqdict(path, pattern=''):
    print 'Generating seqdict...'
    seqdict = {}
    for fn in glob.glob(os.path.join(path, pattern)):
        print 'Processing', fn
        seqid = None
        for line in open(fn):
            if line.startswith('>'):
                seqid = '%s|%s' % (line.strip()[1:], os.path.basename(fn).rstrip('.fas')) 
            else:
                seqdict[seqid] = line.strip()                                
    return seqdict


#
def remove_selfmatches(open_fn, out=sys.stdout):
    
    """
    This removes self-matches from a given raw bowtie output.
    """    
    for line in open_fn:
        line_token = line.strip().split('\t')
        if line_token[0] != line_token[2]:
            out.write('%s\n' % line)    
    pass

"""
def get_all_pairs(open_fn, out=sys.stdout):
    pairs = set([])
    for line in open_fn:
        line_token = line.strip().split('\t')
        pair = tuple(sorted([line_token[0], line_token[2]]))
        pairs.add(pair)
    return list(pairs)

def find_unique_pairs(open_fn, refpairs, out=sys.stdout):
    for line in open_fn:
        line_token = line.strip().split('\t')
        pair = tuple(sorted([line_token[0], line_token[2]]))
        
    pass
"""

""" 
[aligned_read, reference_strand_aligned_to, reference_seq, 5'-match_(0-based), read_sequence]
@HWI-ST377:150:C1116ACXX:5:1101:10126:70489|CTS1_miRNA_R1.unpaired      +       @HWI-ST377:150:C1116ACXX:6:2305:13989:30984|CTS1_snRNA_R1.paired       0       GAGACGTGATGAACCCACTAATTGGTCCGTGTTTCTGATTACCGTGACTGATAAACCTGTGTTTTTCTGATCTTGGAATTCTCGGGTGCCAAGGAACTCC   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII   1473
"""

def get_unique_sequences(open_fn, seqdict, out=sys.stdout):
    seqs = defaultdict(set)
    print 'Getting unique sequences...'
    for line in open_fn:
        #print 
        #print '>>>>>', line
        #print 
        line_token = line.strip().split('\t')
         
        if line_token[0] == line_token[2]:
            # ignore self-matches
            continue
        
        read_seq = line_token[4]
        try:
            ref_seq = seqdict[line_token[2]]
        except:
            'ERROR: seqid', line_token[2], 'not found in seqdict. Aborting.'
            sys.exit(1)            
        
        seqs[ref_seq].add((0, line_token[2]))        
            
        if len(read_seq) ==  len(ref_seq):
            # 100% identical sequences, since data is from 0-error matching
            seqs[ref_seq].add((0, line_token[0]))            
        else:
            # read_seq is contained in ref_seq
            seqs[ref_seq].add((1, line_token[0]))
    
    return seqs
            
        
        
        
        
        
    


def main(argv):
    # remove_selfmatches(open(argv[0]))
    print 'Beeep Beeep'    
    ngspath = '/home/schudoma/ngstemp/cschudoma'
    projectpath = os.path.join(ngspath, 'Sly-Cus_parasitic_RNA-Seq_Kragler')    
    readpath = os.path.join(projectpath, 'fasta')
    workpath = os.path.join(projectpath, 'nodupes')
    workfile = os.path.join(workpath, 'CTS_CTS_noself.map')
        
    seqdict = make_seqdict(readpath, pattern='CTS*.fas')
    
    # test
    # print seqdict['@HWI-ST377:150:C1116ACXX:5:1101:10126:70489|CTS1_miRNA_R1.unpaired']
    out = open(os.path.join(workpath, 'CTS_unique_sequences.fas2'), 'wb') #sys.stdout
    unique_seqs = get_unique_sequences(open(workfile), seqdict)
    
    for seq, hits in unique_seqs.items():
        hits = sorted(list(hits))
        out.write('>%s|%i\n%s\n' % (hits[0][1].split('|')[0], len(hits), seq))
        # break
        
    
    
    
    pass


if __name__ == '__main__': main(sys.argv[1:])

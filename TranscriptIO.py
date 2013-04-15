#!/usr/bin/env python
'''
Created on Apr 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import cPickle as pickle

from Transcript import Transcript
from SNP_Position import SNP_Position

def load_transcript_data(open_fn, sample=''):
    transcript_d = pickle.load(open_fn)
    for k in transcript_d:
        transcript_d[k].calculate_binom_score(sample=sample)
        transcript_d[k].check_for_mobility(sample=sample)
    return transcript_d
        

def read_transcript_data(open_fn):
    transcript_d = {}
    snp_d = {}
    for line in open_fn:
        line = line.strip().split('\t')
        transcript, snp = Transcript(line[:9]), SNP_Position(line[9:])
        
        snp_key = (transcript.contig, snp.pos_1)
        if snp_key in snp_d:
            snp = snp_d[snp_key][0]
        else:
            snp_d[snp_key] = [snp]
        snp_d[snp_key].append(transcript.id_)
        snp.count += 1    
                
        if transcript.id_ not in transcript_d:
            transcript_d[transcript.id_] = transcript
        transcript_d[transcript.id_].snps.append(snp)
        
    return transcript_d, snp_d


def main(argv):
    pass

if __name__ == '__main__': main(sys.argv[1:])

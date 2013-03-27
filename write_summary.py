#!/usr/bin/env python
'''
Created on Mar 27, 2013

@author: Chris
'''

import sys
import os

import pickle

from analyse_snps_by_transcript import Transcript
from analyse_snps_by_transcript import SNP_Position 



def gather_transcript_data(transcripts, transcript_read_data, sample_name):
    for k in transcripts:
        key = transcripts[k].get_key()
        if key not in transcript_read_data:
            transcript_read_data[key] = {}
        
        sample = None
        if 'Col' in sample_name or 'col' in sample_name:
            sample = 'Col'
        elif 'Ped' in sample_name or 'ped' in sample_name:
            sample = 'Ped'
        
        has_covered_snps = reduce(lambda x,y:x or y, [snp.is_covered for snp in transcripts[k].snps])
        if has_covered_snps:
            support_col, total, support_ped, support_both, has_conflict, is_mobile, gmean_pp = transcripts[k].check_snps(sample=sample)
            transcript_read_data[key][sample_name] = total, support_col, support_ped, support_both, has_conflict, is_mobile, gmean_pp
        else:
            transcript_read_data[key][sample_name] = None
    pass

def write_summary(transcript_data, samplenames):
    
    sys.stdout.write('\t'.join(['AGI', 'Contig', 'Type', 'Start', 'End', 'Strand'])) + '\t'
    for sn in samplenames:
        sys.stdout.write('\t'.join(['%s:%s' % (sn, head) 
                                    for head in ['Total_SNPs', 'SNPs_Col', 'SNPs_Ped', 'SNPs_both', 'Conflict?', 'Mobile_Candidate?', 'P']])) + '\t'   
    
    for k in sorted(transcript_data):
        line = '\t'.join(map(str, k)) + '\t'
        for sn in samplenames:
            line += '\t'.join(map(str, transcript_data[k][sn]))
        sys.stdout.write(line + '\n')
    pass
            
        


def main(argv):
    
    transcript_data = {}
    samplenames = [os.path.basename(arg) for arg in argv]
    for fn in argv:
        data = pickle.load(open(fn, 'rb'))
        gather_transcript_data(data, transcript_data, fn)
    
    write_summary(transcript_data, samplenames)
    
    pass



if __name__ == '__main__': main(sys.argv[1:])

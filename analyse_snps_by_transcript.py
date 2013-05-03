#!/usr/bin/env python
'''
Created on Mar 14, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import cPickle as pickle
import time
import datetime

import pysam

from Transcript import Transcript
from SNP_Position import SNP_Position
from TranscriptIO import read_transcript_data, load_transcript_data, get_timestamp 
from count_bases import count_bases


        
def filter_ambiguous_reads(reads):
    all_reads = {}
    # flag all reads whether they are found supporting col (1), ped(2), or cause a conflict
    # because they support both col and ped at different positions (1 | 2 = 3)
    for read in reads:         
        all_reads[read[0]] = all_reads.get(read[0], 0) | (int(read[-1]) + 1)
    # now remove all occurrences of conflicting reads
    # as well as reads that map different SNPs within the same gene (by making sets)
    return [read for read in all_reads if all_reads[read] != 3]
      
def process_pileups(bamfile, snp_d, read_checklist):
    N = len(snp_d)
    i = 1.0    
    for k in snp_d:
        if (i*100/N) % 5 == 0: 
            sys.stderr.write('%s: %i SNPs processed (%i%%)\n' % (get_timestamp(), i, i/N * 100))
        i += 1.0
        snp = snp_d[k][0]
        snp.test += 1
        base_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 
                      'low_qual': 0, 'del': 0, 'mates_disagree': 0, 'unique_reads': 0} 
        reads_at_pos = dict(zip(base_count.keys(), [[] for key in base_count.keys()]))
        valid_reads = []
        # print   snp.get_region(k[0])  
        for col in bamfile.pileup(region=snp.get_region(k[0])):
            if snp.pos_0 == col.pos:
                base_count, weighted, valid_reads, reads_at_pos = count_bases(col)
                break
        if base_count is None or reads_at_pos is None:
            print 'Skipping', snp.get_region(k[0]) 
            continue
        snp.add_pileup(base_count, reads_at_pos)
        
        # print 'V', valid_reads
        valid_reads = set([(read, base==snp.mutation)
                           for read, base, seq, qpos in valid_reads            
                           if base in (snp.refbase, snp.mutation)])
        read_checklist.update(valid_reads)
    
    # This is a simulation to see how ambiguous SNPs can be removed from the read count
    # print snp_d.keys()
    # snp_d[('Chr4', 10952350)][0].r_support_mutation.append('HWI-ST377:127:D0PHGACXX:7:1306:11252:6941')
    # read_checklist.add(('HWI-ST377:127:D0PHGACXX:7:1306:11252:6941', True))
    # print 'x', sorted(snp_d[('Chr4', 10952350)][0].r_support_refbase[:10])
    # print 'x', sorted(snp_d[('Chr4', 10951240)][0].r_support_refbase[:10])
    
    filtered_reads = set(filter_ambiguous_reads(read_checklist))
    # print len(read_checklist), len(filtered_reads) 
    for k in snp_d:
        snp = snp_d[k][0]
        if snp.is_covered:
            snp.cleanse_counts(filtered_reads)
        
    pass

def write_data(transcript_d, sample='', prefix=''):
    fo = open(prefix + 'SNPTABLE.csv', 'w')
    candidates = []
    sys.stderr.write('%s: Writing SNP/transcript table...\n' % (get_timestamp())) 
    for k in sorted(transcript_d):
        transcript_shown = False
        # show_snps = reduce(lambda x,y:x or y, [snp.is_covered for snp in transcript_d[k].snps])
        if not transcript_d[k].has_covered_snps(): 
            continue         
        # transcript_d[k].calculate_binom_score(sample=sample)
        # transcript_d[k].check_for_mobility(sample=sample)
        if transcript_d[k].is_mobile and transcript_d[k].binom_score > 10.0 ** -5:
            candidates.append(k)       
        for snp in transcript_d[k].snps:
            if not snp.is_covered:
                continue
            if not transcript_shown:
                fo.write(transcript_d[k].get_full_string() + '\n')
                transcript_shown = True
            fo.write('\t' + str(snp) + '\n')
    fo.close()
    fo = open(prefix + 'RANKED.csv', 'w')
    sys.stderr.write('%s: Ranking candidate mobile transcripts (%i candidates)...\n' % (get_timestamp(), len(candidates)))
    ranked = sorted([transcript_d[k] for k in candidates if transcript_d[k].is_mobile], 
                    key=lambda x:x.binom_score, reverse=True)
    sys.stderr.write('%s: Writing ranked transcripts...\n' % (get_timestamp()))
    for transcript in ranked:
        fo.write(transcript.get_full_string() + '\n')
    fo.close()
    
    sys.stderr.write('%s: Generating candidate shortlist...\n' % get_timestamp())
    snp_ratio_cutoff = 0.5
    if sample.lower() == 'col':    
        shortlist = [transcript for transcript in ranked 
                     if float(transcript.support_ped)/transcript.support_total >= snp_ratio_cutoff]
    elif sample.lower() == 'ped':
        shortlist = [transcript for transcript in ranked 
                     if float(transcript.support_col)/transcript.support_total >= snp_ratio_cutoff]
    n_snp_cutoff = 3
    shortlist = [transcript for transcript in shortlist
                 if transcript.support_total >= n_snp_cutoff]
    fo = open(prefix + 'SHORTLIST.csv', 'w')
    sys.stderr.write('%s: Writing candidate shortlist...\n' % get_timestamp())
    for transcript in shortlist:
        fo.write(transcript.get_full_string() + '\n')
    fo.close()
    pass
    
    

def show_data(transcript_d, snp_d, sample=None):
    for k in sorted(transcript_d):
        transcript_shown = False
        # print transcript_d[k]
        show_snps = reduce(lambda x,y:x or y, [snp.is_covered for snp in transcript_d[k].snps])
        if not show_snps:
            continue        
        transcript_d[k].calculate_binom_score(sample=sample)
        transcript_d[k].check_for_mobility(sample=sample)
        for snp in transcript_d[k].snps:
            if not snp.is_covered:
                continue
            if not transcript_shown:
                print transcript_d[k].get_full_string()
                transcript_shown = True
            print '\t' + str(snp)
            
    pass    



def main(argv):    
    """
    Input EITHER (to generate per-transcript-results from pileup counts)
    0. transcript/snp data
    1. sample name {Col,Ped}
    2. output prefix (output is <prefix>_<DATATYPE>.pickled)
    3-n. number of bamfiles
    OR (to summarise per-transcript-results)
    0. name of a <arbitrary>_TRANSCRIPTDATA.pickled file 
    1. sample name {Col,Ped}
    """
    
    if len(argv) > 3:
        sys.stderr.write('%s: MODE1\n' % get_timestamp())
        sys.stderr.write('')
        transcript_d, snp_d = read_transcript_data(open(argv[0]))
        # show_data(transcript_d, snp_d)
        sample = argv[1]
        prefix = argv[2]
        read_checklist = set([])
        for bam_fn in argv[3:]:  
            sys.stderr.write('%s: Processing file %s...\n' % (get_timestamp(), bam_fn))      
            process_pileups(pysam.Samfile(bam_fn, 'rb'), snp_d, read_checklist)    
            # show_data(transcript_d, snp_d)
            # break
        # print list(read_checklist), len(read_checklist)
        ts = get_timestamp()
        pickle.dump(snp_d, open(prefix + '_SNPDATA.pickled', 'wb'))
        pickle.dump(transcript_d, open(prefix + '_TRANSCRIPTDATA.pickled', 'wb'))
        pickle.dump(read_checklist, open(prefix + '_READCHECKLIST.pickled', 'wb'))
	# sys.exit(0)
    else:
        sys.stderr.write('%s: MODE2\n' % get_timestamp())
        sys.stderr.write('%s: Loading data from %s.\n' % (get_timestamp(), argv[0]))
        # transcript_d = pickle.load(open(argv[0], 'rb'))
        transcript_d = load_transcript_data(open(argv[0], 'rb'), sample=argv[1])
        sys.stderr.write('%s: Finished loading %s.\n' % (get_timestamp(), argv[0]))
        prefix = argv[0].rstrip('TRANSCRIPTDATA.pickled')
    
    # print '======'
    
    #show_data(transcript_d, None, sample=argv[1])
    write_data(transcript_d, sample=argv[1], prefix=prefix)
    pass

if __name__ == '__main__': main(sys.argv[1:])

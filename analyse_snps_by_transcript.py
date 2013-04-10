#!/usr/bin/env python
'''
Created on Mar 14, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import cPickle as pickle

import math
import time
import datetime

import scipy.stats as stats
import scipy

from collections import Counter

import pysam

TRANSCRIPT_ATTR = ['contig', None, 'type', 'start', 'end', None, 'strand', None, 'id_']
TRANSCRIPT_CASTS = [str, None, str, int, int, None, str, None, str]

PILEUP_ATTR = ['support_refbase', 'support_mutation', 'N_at_pos', 'mismatches', 'lowqual', 'indel', 'n_aligned_reads']
R_PILEUP_ATTER = ['r_support_refbase', 'r_support_mutation', 'r_N_at_pos', 'r_mismatches', 'r_lowqual', 'r_indel', 'r_n_aligned_reads']

from count_bases import count_bases

def stats_binom_wrapper(x, n, p):
    # 0.13.0.dev-2bd1af0
    if scipy.__version__ >= '0.13' or (p != 1 and p != 0):
        return stats.binom.pmf(x, n, p)
    else:
        return float(p)


class Transcript(object):
    def __init__(self, attributes):
        self.snps = []
        for name, val, cast in [x 
                                for x in zip(TRANSCRIPT_ATTR, attributes, TRANSCRIPT_CASTS) 
                                if not x[0] is None]:
            if name == 'id_':
                val = val.split(';')[0].lstrip('ID=')
                pass
            setattr(self, name, cast(val))
        pass
    def __str__(self):
        return '\t'.join(map(str, [self.id_, self.contig, self.type, self.start, self.end, '"%s"' % self.strand]))
    
    def get_key(self):
        return (self.id_, self.contig, self.type, self.start, self.end, self.strand)
    
    def get_full_string(self):
        string = str(self) + '\t\tCOL: %i/%i, PED: %i/%i, BOTH: %i/%i, CONFLICT=%s, MOBILE_CANDIDATE=%s, P=%s'
        return string % (self.support_col, self.support_total, self.support_ped, self.support_total,
                         self.support_both, self.support_total,
                         self.has_conflict, self.is_mobile, str(self.binom_score))
    
    def calculate_binom_score(self, sample=''):
        col_total, ped_total = 1.0, 1.0
        
        for snp in self.snps:
            if snp.is_covered:                
                col_total += len(snp.r_support_refbase)
                ped_total += len(snp.r_support_mutation)
        both_total = col_total + ped_total
        
        all_p_pileup = []
        for snp in self.snps:
            if snp.is_covered:
                n_col, n_ped = float(len(snp.r_support_refbase) + 1), float(len(snp.r_support_mutation) + 1)
                if both_total == 0:
                    snp.p_pileup = 666.0
                    continue
                if sample.lower() == 'col':
                    p_pileup = stats_binom_wrapper(n_col, n_col + n_ped, col_total/both_total)
                elif sample.lower() == 'ped':
                    p_pileup = stats_binom_wrapper(n_ped, n_col + n_ped, ped_total/both_total)
                else:
                    sys.stderr.write('WRONG SAMPLE:%s. Exiting.\n' % str(sample))
                    sys.exit(1)
                snp.p_pileup = p_pileup
                all_p_pileup.append(p_pileup)        
        
        self.binom_score = stats.mstats.gmean(all_p_pileup)
        pass
    

    def check_for_mobility(self, sample='', sampledict={'col': 'refbase', 'ped': 'mutation'}):
        numerator = sampledict.get(sample.lower(), 'refbase')        
        self.has_conflict = False
         
        for snp in self.snps:
            if snp.is_covered:
                if snp.get_read_ratio(numerator=numerator) >= 0.5:                
                    snp.ratio_ok = True
                else:
                    snp.ratio_ok = False
                    self.has_conflict = True
        
        c = Counter([snp.check_support() for snp in self.snps if snp.is_covered])
        self.support_col, self.support_ped = c['Col'] + c['both'], c['Ped'] + c['both']
        self.support_total = sum(c.values())
        self.support_both = c['both']
        
        if self.has_conflict:
            self.is_mobile = False        
        elif sample.lower() == 'col':
            self.is_mobile = self.support_ped > 0
        elif sample.lower() == 'ped':
            self.is_mobile = self.support_col > 0
        else:
            self.is_mobile = False
                        
        pass
    
    pass
    
            
class SNP_Position(object):
    def __init__(self, attributes):        
        self.pos_0 = int(attributes[0])
        self.pos_1 = int(attributes[1])
        self.count = 0
        self.test = 0
        self.is_covered = False
        for name, val in [x.split('=') for x in attributes[2].split(';')]:
            try:
                val = float(val)
            except:
                pass
            setattr(self, name, val)
        pass
    
    def __str__(self):
        
        if hasattr(self, 'support_refbase'):            
            
            STR_RATIO_OK = ''
            if hasattr(self, 'ratio_ok') and not self.ratio_ok:
                STR_RATIO_OK = '!!'
            STR_PPILEUP = ''
            if hasattr(self, 'p_pileup'):
                STR_PPILEUP =  str(self.p_pileup)
            readstats = '\t'.join(map(str, ['', '', #self.support_refbase, self.support_mutation,
                                            len(self.r_support_refbase), len(self.r_support_mutation),                                            
                                            'RELIABLE_SNP=%s' % (self.RELIABLE == 'YES'),
                                            self.check_support(),
                                            '%.5f' % self.get_read_ratio() + ('%s' % STR_RATIO_OK),
                                            STR_PPILEUP 
                                            ]))
        else:
            readstats = '\t'.join((['0'] * 8) + ['N/A']) #'0\t0\t0\t0\tN/A'        
        return '\t'.join(map(str, [self.pos_0, self.pos_1, self.refbase, self.mutation, readstats]))
    
    
    
    
    
    
    def get_region(self, contig):
        return '%s:%i-%i' % (contig, self.pos_0, self.pos_1)
    def get_read_ratio(self, numerator='refbase'):
        if numerator == 'refbase':
            return (len(self.r_support_refbase) + 1.0) / (len(self.r_support_refbase) + len(self.r_support_mutation) + 2.0)
        elif numerator == 'mutation':
             return (len(self.r_support_mutation) + 1.0) / (len(self.r_support_refbase) + len(self.r_support_mutation) + 2.0)
        else:
            return 0.0  
    
    def add_pileup(self, pileup_count, reads_at_pos):
        
        for attr in PILEUP_ATTR:
            if not hasattr(self, attr):
                setattr(self, attr, 0)
        for attr in R_PILEUP_ATTER:
            if not hasattr(self, attr):
                setattr(self, attr, [])
        
        self.support_refbase += pileup_count[self.refbase]
        self.support_mutation += pileup_count[self.mutation]
        self.N_at_pos += pileup_count['N']
        self.mismatches += sum([pileup_count[c] for c in 'ACGT']) - \
        (pileup_count[self.refbase] + pileup_count[self.mutation])
        self.lowqual += pileup_count['low_qual']
        self.indel += pileup_count['del']
        self.n_aligned_reads += pileup_count['unique_reads'] - self.indel
        
        self.r_support_refbase += reads_at_pos[self.refbase]
        self.r_support_mutation += reads_at_pos[self.mutation]
        self.r_N_at_pos += reads_at_pos['N']
        
        if len(self.r_support_refbase) > 0 or len(self.r_support_mutation) > 0:
            self.is_covered = True
        
        abc = list(set(list('ACGT')) - set([self.refbase, self.mutation]))
        for c in abc:
            self.r_mismatches += reads_at_pos[c]
        self.r_lowqual += reads_at_pos['low_qual']
        self.r_indel += reads_at_pos['del']
        self.r_n_aligned_reads += list(set(reads_at_pos['unique_reads']) - set(self.r_indel))        
        pass
    def check_support(self, cutoff=3, trd={0: 'N/A', 1: 'Ped', 2: 'Col', 3: 'both'}):
        try:
            flag = (int(len(self.r_support_refbase) >= cutoff) << 1) | int(len(self.r_support_mutation) >= cutoff)
        except:
            # print self.pos_0, self.pos_1
            # print dir(self)
            flag = 0
        if trd is not None:
            return trd[flag]
        else: 
            return flag 
        
    def get_read_ratio2(self, cutoff=3, sample='refbase'):
        n_refbase = len(self.r_support_refbase)
        if n_refbase < cutoff: 
            n_refbase = 0
        n_mutation =  len(self.r_support_mutation)
        if n_mutation < cutoff:
            n_mutation = 0
        
        if sample == 'refbase':
            ratio = float(n_refbase + 1) / float(n_mutation + 1)
        else:
            ratio =  float(n_mutation + 1) / float(n_refbase + 1)
        
        return ratio, ratio/(n_refbase + n_mutation + 2)    
        
    
    def cleanse_counts(self, valid_reads):        
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) - set(self.r_support_refbase))
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) - set(self.r_support_mutation))                
        self.r_support_refbase = list(set(self.r_support_refbase).intersection(valid_reads))
        self.r_support_mutation = list(set(self.r_support_mutation).intersection(valid_reads))        
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) | set(self.r_support_refbase))
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) | set(self.r_support_mutation))
        pass
                
    pass
      

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
        show_snps = reduce(lambda x,y:x or y, [snp.is_covered for snp in transcript_d[k].snps])
        if not show_snps:
            continue         
        transcript_d[k].calculate_binom_score(sample=sample)
        transcript_d[k].check_for_mobility(sample=sample)
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

def get_timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')


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
        transcript_d = pickle.load(open(argv[0], 'rb'))
        sys.stderr.write('%s: Finished loading %s.\n' % (get_timestamp(), argv[0]))
        
    
    # print '======'
    
    #show_data(transcript_d, None, sample=argv[1])
    write_data(transcript_d, sample=argv[1], prefix=argv[0].rstrip('TRANSCRIPTDATA.pickled'))
    pass

if __name__ == '__main__': main(sys.argv[1:])

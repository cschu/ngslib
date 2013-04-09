#!/usr/bin/env python
'''
Created on Mar 14, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import pickle

import math
import time
import datetime

import scipy.stats as stats

from collections import Counter

import pysam

TRANSCRIPT_ATTR = ['contig', None, 'type', 'start', 'end', None, 'strand', None, 'id_']
TRANSCRIPT_CASTS = [str, None, str, int, int, None, str, None, str]

PILEUP_ATTR = ['support_refbase', 'support_mutation', 'N_at_pos', 'mismatches', 'lowqual', 'indel', 'n_aligned_reads']
R_PILEUP_ATTER = ['r_support_refbase', 'r_support_mutation', 'r_N_at_pos', 'r_mismatches', 'r_lowqual', 'r_indel', 'r_n_aligned_reads']

from count_bases import count_bases

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
    

    def check_snps(self, cutoff=3, sample=None):
        
        def stats_binom_wrapper(x, n, p):
            if p == 1: # and x == n:
                return 1.0
            elif p == 0: # and x == 0:
                return 0.0
            else:
                return stats.binom.pmf(x, n, p)
            
        
        read_ratios = []
        
        if sample == 'Col':
            numerator = 'refbase'
        else:
            numerator = 'mutation'
            
        col_total, ped_total = 0.0, 0.0
        for snp in self.snps:
            if snp.is_covered:
                read_ratios.append(snp.get_read_ratio(numerator=numerator))
                col_total += len(snp.r_support_refbase)
                ped_total += len(snp.r_support_mutation)
        both_total = col_total + ped_total
        
        # print 'C=', col_total, 'P=', ped_total, 'B=', both_total
             
        """
        cdf(x, n, p, loc=0)    Cumulative density function.
        x : array_like quantiles
        n, p : array_like shape parameters
        loc : loc : array_like, optional location parameter (default=0)
        """
        all_p_pileup = []
        for snp in self.snps:
            if snp.is_covered:
                n_col, n_ped = float(len(snp.r_support_refbase)), float(len(snp.r_support_mutation))
                if both_total == 0:
                    snp.p_pileup = 666.0
                    continue
                if sample == 'Col':
                    # print '**'
                    # print n_col, n_ped, col_total, both_total, col_total/both_total
                    # p_pileup = stats.binom.pmf(n_col + 1, n_col + n_ped + 2, col_total/both_total)
                    p_pileup = stats_binom_wrapper(n_col + 1, n_col + n_ped + 2, col_total/both_total)
                    # print '***'
                elif sample == 'Ped':
                    # p_pileup = stats.binom.pmf(n_ped, n_col + n_ped, ped_total/both_total)
                    p_pileup = stats_binom_wrapper(n_ped, n_col + n_ped, ped_total/both_total)
                else:
                    sys.stderr.write('WRONG SAMPLE:%s. Exiting.\n' % str(sample))
                    sys.exit(1)
                # if p_pileup >= 0.5:
                #    p_pileup = 1.0 - p_pileup
                snp.p_pileup = p_pileup
                all_p_pileup.append(p_pileup)
        
        
        # print all_p_pileup
        gmean_p_pileup = stats.mstats.gmean(all_p_pileup)
        
                
        
        has_conflict = False
        #if len(read_ratios) >= 0:
            # mean = sum(read_ratios) / float(len(read_ratios))
            # stdev = math.sqrt(sum([(ratio - mean) ** 2 for ratio in read_ratios]) / float(len(read_ratios) - 1))
        
         
        for snp in self.snps:
            if snp.is_covered:
                if snp.get_read_ratio(numerator=numerator) >= 0.5:
                #if (mean - stdev) <= snp.get_read_ratio(numerator=numerator) <= (mean + stdev):
                    snp.ratio_ok = True
                else:
                    snp.ratio_ok = False
                    has_conflict = True
        
        c = Counter([snp.check_support() for snp in self.snps if snp.is_covered])
        support_col, support_ped = c['Col'] + c['both'], c['Ped'] + c['both']
        total = sum(c.values())
        
        if sample == 'Col':
            is_mobile = support_ped > 0
        elif sample == 'Ped':
            is_mobile = support_col > 0
        else:
            is_mobile = False
        
        string = '\tCOL: %i/%i, PED: %i/%i, BOTH: %i/%i' % (support_col, total, support_ped, total, c['both'], total)
        string += ', CONFLICT=%s' % has_conflict
        string += ', MOBILE_CANDIDATE=%s' % is_mobile        
        string += ', P=%s' % str(gmean_p_pileup)
        
        #return string
        return support_col, total, support_ped, c['both'], has_conflict, is_mobile, gmean_p_pileup
    
    def check_snps_old(self, cutoff=3, sample=None):
        support_flags = [snp.check_support() for snp in self.snps if snp.is_covered]
        c = Counter(support_flags)
        no_support = c['-']
        del c['-']
        support_col = c['Col'] + c['both']
        support_ped = c['Ped'] + c['both']
        support_both = c['both']
        total = sum(c.values()) + no_support
        # del c['both']
        # total = sum(c.values())
        
        string = 'COL: %i/%i, PED: %i/%i, BOTH: %i/%i' % (support_col, total, support_ped, total, c['both'], total)
        # mc = c.most_common(1)[0]
        #if mc[0] == 'Col' and mc[1]
        
        is_mobile, has_conflict = False, False
        if sample == 'Col':
            is_mobile = support_ped > 0
        elif sample == 'Ped':
            is_mobile = support_col > 0
        
        has_conflict = (support_ped > 0 and support_col > 0 and support_ped != support_col)
        
        # string += ', CONFLICT=%s' % (c['Col'] > 0 and c['Ped'] > 0)
        string += ', CONFLICT=%s' % has_conflict
        
        # is_mobile = (sample == 'Col' and support_ped > 0) or (sample == 'Ped' and support_col > 0)
        string += ', MOBILE=%s' % is_mobile
        
        string += '\t' + str(c)
        
        
        return string
        
        
                
        
    
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
            readstats = '\t'.join(map(str, ['', '', #self.support_refbase, self.support_mutation,
                                            len(self.r_support_refbase), len(self.r_support_mutation),                                            
                                            'RELIABLE_SNP=%s' % (self.RELIABLE == 'YES')]))
        else:
            readstats = '0\t0\t0\t0\tN/A'        
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



"""
string = '\tCOL: %i/%i, PED: %i/%i, BOTH: %i/%i' % (support_col, total, support_ped, total, c['both'], total)
        string += ', CONFLICT=%s' % has_conflict
        string += ', MOBILE_CANDIDATE=%s' % is_mobile        
        string += ', P=%s' % str(gmean_p_pileup)
        
        #return string
        return support_col, total, support_ped, c['both'], has_conflict, is_mobile, gmean_p_pileup
"""


def show_data(transcript_d, snp_d, sample=None):
    for k in sorted(transcript_d):
        transcript_shown = False
        # print transcript_d[k]
        show_snps = reduce(lambda x,y:x or y, [snp.is_covered for snp in transcript_d[k].snps])
        if not show_snps:
            continue        
        for snp in transcript_d[k].snps:
            if not snp.is_covered:
                continue
            if not transcript_shown:
                support_col, total, support_ped, support_both, has_conflict, is_mobile, gmean_pp = transcript_d[k].check_snps(sample=sample)
                string = '\t\tCOL: %i/%i, PED: %i/%i, BOTH: %i/%i, CONFLICT=%s, MOBILE_CANDIDATE=%s, P=%s'                
                print transcript_d[k], string % (support_col, total, support_ped, total, support_both, total, has_conflict, is_mobile, str(gmean_pp))
                transcript_shown = True
            print '\t',
            # print snp, '\t', 'COVERED=%s' % snp.is_covered, '\t', snp.check_support(),
            print snp, '\t', snp.check_support(),
            print '\t', '%.5f' % snp.get_read_ratio(),
            if hasattr(snp, 'ratio_ok') and not snp.ratio_ok:
                print '!!', 
            else:
                print '',
            if hasattr(snp, 'p_pileup'):
                print '\t', str(snp.p_pileup)
            else:
                print ''
            # print
#    for k in sorted(snp_d):
#        print k
#        for tx in snp_d[k][1:]:
#            print '\t' + tx
    pass    

def get_timestamp():
    return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')


def main(argv):
    
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
    
    show_data(transcript_d, None, sample=argv[1])
    pass

if __name__ == '__main__': main(sys.argv[1:])

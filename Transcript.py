#!/usr/bin/env python
'''
Created on Apr 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
from collections import Counter

import scipy.stats as stats
import scipy

TRANSCRIPT_ATTR = ['contig', None, 'type', 'start', 'end', None, 'strand', None, 'id_']
TRANSCRIPT_CASTS = [str, None, str, int, int, None, str, None, str]


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
        string = '%s\t%s\t%s:%i-%i:%s' % (self.id_, self.type, self.contig, self.start, self.end, self.strand)
        string += '\tSNPs=%i' % len(self.snps)
        return string
        # return '\t'.join(map(str, [self.id_, self.contig, self.type, self.start, self.end, '"%s"' % self.strand]))
    
    def get_key(self):
        return (self.id_, self.contig, self.type, self.start, self.end, self.strand)
    
    def has_covered_snps(self):
        return reduce(lambda x,y:x or y, [snp.is_covered for snp in self.snps])
    
    def get_full_string(self):
        string = str(self) + '\t\tCOL: %i/%i, PED: %i/%i, BOTH: %i/%i, CONFLICT=%s, MOBILE_CANDIDATE=%s, P=%s'
        # total = self.support_total
        total = len(self.snps)        
        return string % (self.support_col, total, self.support_ped, total,
                         self.support_both, total,
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

        # this applies Default Filter Rule 2         
        for snp in self.snps:
            if snp.is_covered:
                if snp.get_read_ratio(numerator=numerator) >= 0.5:                
                    snp.ratio_ok = True
                else:
                    snp.ratio_ok = False
                    self.has_conflict = True
        
        # this applies Default Filter Rule 1 
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


def main(argv):
    pass

if __name__ == '__main__': main(sys.argv[1:])

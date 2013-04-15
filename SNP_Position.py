#!/usr/bin/env python
'''
Created on Apr 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

PILEUP_ATTR = ['support_refbase', 'support_mutation', 'N_at_pos', 'mismatches', 'lowqual', 'indel', 'n_aligned_reads']
R_PILEUP_ATTER = ['r_support_refbase', 'r_support_mutation', 'r_N_at_pos', 'r_mismatches', 'r_lowqual', 'r_indel', 'r_n_aligned_reads']

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



def main(argv):
    pass

if __name__ == '__main__': main(sys.argv[1:])

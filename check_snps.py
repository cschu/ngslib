#!/usr/bin/python

import os
import re
import sys
import math

from snptools import SNPLine, SNPDict

###
def do_something():
    return None 


def read_snp_data(open_fn):
    snp_d = {}
    for line in open_fn:
        snpline = SNPLine(line, lineid=len(snp_d) + 1)
        snp_d[(snpline.seqid, snpline.pos)] = snpline
    return snp_d


SNP_CASTS = {
    'contig': str,
    'pos': int,
    'ref_base': str,
    'A_count': int,
    'T_count': int,
    'G_count': int,
    'C_count': int,
    'bad_count': int,
    'total_count': int
    }

def dic2str(dic, sep='/'):
    return sep.join(['%s:%s' % tuple(map(str, item))
                     for item in dic.items()])

class SNPPileup(object):
    def __init__(self, fields, values, casts):
        self.basecounts = {}
        for field, value in zip(fields, values):
            if not re.match('[ACGT]_count', field):
                setattr(self, field, casts[field](value))
            else:
                self.basecounts[field[0]] = casts[field](value)                
        pass    
    def __repr__(self):
        

        outstr = [self.contig, self.pos, self.ref_base, 
                  dic2str(self.basecounts),
                  sum(self.basecounts.values()) + self.bad_count]
        return '\t'.join(map(str, outstr))
    
    pass

def read_snp_pileups(open_fn, casts=SNP_CASTS):
    pileups_d = {}
    fields = open_fn.readline().strip().split('\t')
    for line in open_fn:
        snpline = SNPPileup(fields, line.split('\t'), casts)
        pileups_d[(snpline.contig, snpline.pos)] = snpline
    return pileups_d


"""
data:
Ped-0   1       125155  A       G       38      9       1.0     1.0
Ped-0   1       125158  C       A       40      9       1.0     1.0
Ped-0   1       125344  T       G       40      17      1.0     1.0

pileups:
contig  pos     ref     A_count T_count G_count C_count bad     total
Chr1    125155  A       0       0       2       0       0       2
Chr1    125158  C       2       0       0       0       0       2
Chr1    186861  A       0       0       1       0       0       1
"""

###
def main(argv):

    # snp_data = # read_snp_data(open(argv[0]))
    snp_data = SNPDict(open(argv[0]))
    snp_pileups = read_snp_pileups(open(argv[1]))

    min_reads = 10
    min_reads_with_snp = 5

    base_ecotype = 'col-0'
    snp_ecotype = 'ped-0'    

    for snp_id, snp in sorted(snp_data.items()):
        if snp.refbase not in 'ACGT': 
            continue
        if snp.mutation not in 'ACGT':
            continue
        if snp_id not in snp_pileups:
            continue
        if snp_pileups[snp_id].total_count < min_reads:
            continue
        if snp.mutation == '-': 
            continue
        if snp_pileups[snp_id].basecounts[snp.mutation] < min_reads_with_snp:
            continue

        sorted_counts = sorted(snp_pileups[snp_id].basecounts.items(),
                               key=lambda x:x[1])
        mutation_count = snp_pileups[snp_id].basecounts[snp.mutation]
        ref_count = snp_pileups[snp_id].basecounts[snp.refbase]
        """
        if sorted_counts[-1][1] == mutation_count:
            if mutation_count >= 1.5 * sorted_counts[-2][1]:
                # SNP is highest count
                comment = snp_ecotype
            else:
                comment = 'unknown'
        else:
            comment = base_ecotype
        """
        comment = ''
        if sum([x[1] for x in sorted_counts[:2]]) > 0:
            comment = 'multi_snp'
        else:
            valid_counts = tuple(sorted([mutation_count, ref_count]))
            max_counts = tuple([x[1] for x in sorted_counts[-2:]])
            # print valid_counts, max_counts, valid_counts == max_counts
            if valid_counts == max_counts:
                # print float(max_counts[-2])/max_counts[-1]
                if float(max_counts[-2])/max_counts[-1] <= 0.25:
                    if sorted_counts[-1][1] == ref_count:
                        comment = base_ecotype
                    elif sorted_counts[-1][1] == mutation_count:
                        comment = snp_ecotype
                    else:
                        comment = 'unknown'
                else:
                    comment = '%s/%s' % (base_ecotype, snp_ecotype)
                pass
            else:
                comment = 'weird_snp'
            pass

        # Peol-shoot_TGACCA_L007_R1_001.paired.fq.gsnap.paired_mult.bam.samsnp
        read_src = os.path.basename(argv[1]).rstrip('.bam.samsnp')
        # Peol-shoot_TGACCA_L007_R1_001.paired.fq.gsnap.paired_mult
        read_src = re.sub('_[ACGT]{6}_L007_', '_', read_src)
        # Peol-shoot_R1_001.paired.fq.gsnap.paired_mult
        read_src = read_src.replace('.fq.gsnap', '')

        print '\t'.join(map(str,[read_src, 
                                 snp_pileups[snp_id], 
                                 snp, comment]))
        
        
        

            
            



    return None

if __name__ == '__main__': main(sys.argv[1:])

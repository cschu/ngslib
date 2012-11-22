#!/usr/bin/env python
'''
Created on Nov 16, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re

import pysam


def count_bases(col, cutoff=-5):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 'low_qual': 0, 'del': 0}    
    bad_reads = []
    
    for read in col.pileups:
        if read.is_del == 0:
            base = read.alignment.seq[read.qpos]
            qual = ord(read.alignment.qual[read.qpos]) - 33
            if qual > cutoff and base != 'N':
                counts[base] += 1
            elif qual > cutoff:            
                bad_reads.append(('N', read.alignment.qname))
                counts['N'] += 1
            else:
                bad_reads.append(('low_qual', read.alignment.qname))
                counts['low_qual'] += 1
        else:
            bad_reads.append(('del', read.alignment.qname))
            counts['del'] += 1
    # counts['bad'] = len(bad_reads)
    counts['T'] += counts['U']
    del counts['U']
    return counts #, bad_reads
            
            





"""
file: ped-N.unpaired.sorted_hits_snp_intra.gff
Chr1    .       snp     10000334        10000335        .       .       .       
quality=24;supporting_reads=4;concordance=1.0;avg_align_overl_reads=1.0;
refbase=T;mutation=A;gene_start=9996824;gene_end=10000447;gene_strand=+;gene_ID=AT1G28440
"""



def main(argv):
    
    snp_fn = argv[0]
    bam_fn = snp_fn[:snp_fn.find('_')] + '.bam'
    bamfile = pysam.Samfile(bam_fn, 'rb')
    
    gene_d = {}
    
    
    header = ['AGI', 'Pos_SNP', '#Reads', '#Reads_Col', '#Reads_Ped', 
              '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']
    sys.stdout.write(','.join(map(str, header)) + '\n')
    
    for line in open(snp_fn, 'rb'):
        snp = line.strip().split('\t')
        attributes = dict([attr.split('=') for attr in snp[8].split(';')])
        if attributes['refbase'] not in 'ACGTU':
            continue
        if attributes['mutation'] not in 'ACGTU':
            continue
        
        contig, start, end = snp[0], int(snp[3]), int(snp[4])
        region = '%s:%i-%i' % (contig, start, end)
        # print region 
        
        base_count = None
        for col in bamfile.pileup(region=region):
            if start == col.pos:
                base_count = count_bases(col)
        
        if base_count is not None:
            out = [attributes['gene_ID'], end, 
                   sum(base_count.values()),
                   base_count[attributes['refbase']],
                   base_count[attributes['mutation']],
                   base_count['N'],
                   base_count['low_qual'],
                   base_count['del'],
                   sum([base_count[v] for v in 'ACGT']) - \
                   (base_count[attributes['refbase']] + base_count[attributes['mutation']])]
            
            sys.stdout.write(','.join(map(str, out)) + '\n')
                
                
            pass
        
        pass
    bamfile.close()
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python
'''
Created on Nov 16, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re

import pysam

def count_bases(col, cutoff=-5, mult_counts=None):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 'low_qual': 0, 'del': 0}    
    bad_reads = []
    unique_reads = set([])
    
    for read in col.pileups:
        qname = read.alignment.qname
        unique_reads.add(qname)
        try:
            divisor = mult_counts[qname]
        except:
            divisor = 1.0
        increment = 1 / divisor
        
        if read.is_del == 0:    
            base = read.alignment.seq[read.qpos]
            qual = ord(read.alignment.qual[read.qpos]) - 33
            if qual > cutoff and base != 'N':
                counts[base] += increment
            elif qual > cutoff:            
                bad_reads.append(('N', qname))
                counts['N'] += increment
            else:
                bad_reads.append(('low_qual', qname))
                counts['low_qual'] += increment
        else:
            bad_reads.append(('del', qname))
            counts['del'] += increment
    # counts['bad'] = len(bad_reads)
    counts['T'] += counts['U']
    del counts['U']
    counts['unique_reads'] = len(unique_reads)
    return counts #, bad_reads



def analyse_snps(bamfile, snpfile_open, hit_mode, mult_counts=None, out=sys.stdout):
    header = ['', 'Pos_SNP', '#Reads', '#Reads_Col', '#Reads_Ped', 
              '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']
    if hit_mode == 'intra':
        header[0] = 'AGI'
    elif hit_mode == 'inter':
        header[0] = 'Contig'
    else:
        sys.stderr.write('Wrong mode: %s. Exiting.\n' % hit_mode)
        sys.exit(1)
    out.write(','.join(map(str, header)) + '\n')
    
    for line in snpfile_open:
        snp = line.strip().split('\t')
        attributes = dict([attr.split('=') for attr in snp[8].split(';')])
        if attributes['refbase'] not in 'ACGTU':
            continue
        if attributes['mutation'] not in 'ACGTU':
            continue
        contig, start, end = snp[0], int(snp[3]), int(snp[4])
        region = '%s:%i-%i' % (contig, start, end) 
        
        base_count = None
        for col in bamfile.pileup(region=region):
            if start == col.pos:
                base_count = count_bases(col, mult_counts=mult_counts) 
    
        if base_count is not None:            
            out = [None, end, 
                   sum(base_count.values()) - base_count['unique_reads'],
                   # base_count['unique_reads'],
                   base_count[attributes['refbase']],
                   base_count[attributes['mutation']],
                   base_count['N'],
                   base_count['low_qual'],
                   base_count['del'],
                   sum([base_count[v] for v in 'ACGT']) - \
                   (base_count[attributes['refbase']] + base_count[attributes['mutation']])]
            if hit_mode == 'intra':
                out[0] = attributes['gene_ID']
            elif hit_mode == 'inter':
                out[0] = contig
            else:
                sys.stderr.write('Wrong mode: %s. Exiting.\n' % hit_mode)
                sys.exit(1)
            
            sys.stdout.write(','.join(map(str, out)) + '\n')                
            pass
        pass
    bamfile.close()
    pass

    


def main(argv):
    
    snp_fn = argv[0]
    hit_mode = argv[1]
    multi_hits = 'mult' in argv[2:]
    
    bam_fn = snp_fn[:snp_fn.find('_')] + '.bam'
    mult_counts = None
    if multi_hits:
        items = [line.strip().split('\t') 
                 for line in open(bam_fn + '.mult_counts', 'rb').readlines()]
        mult_counts = dict([(item[0], float(item[1])) for item in items])
        pass
    
    analyse_snps(pysam.Samfile(bam_fn, 'rb'), open(snp_fn, 'rb'), hit_mode, mult_counts=mult_counts, out=sys.stdout)
    pass

if __name__ == '__main__': main(sys.argv[1:])
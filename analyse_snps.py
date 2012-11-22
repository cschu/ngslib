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



def analyse_snps(bamfile, snpfile_open, hit_mode, out=sys.stdout):
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
                base_count = count_bases(col) 
    
        if base_count is not None:            
            out = [None, end, 
                   sum(base_count.values()),
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
    bam_fn = snp_fn[:snp_fn.find('_')] + '.bam'
    
    analyse_snps(pysam.Samfile(bam_fn, 'rb'), open(snp_fn, 'rb'), hit_mode, out=sys.stdout)
    pass

if __name__ == '__main__': main(sys.argv[1:])
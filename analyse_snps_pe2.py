#!/usr/bin/env python
'''
Created on Nov 16, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os
import re
import pickle

import pysam


def count_bases(col, cutoff=-5, mult_counts=None):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 
              'low_qual': 0, 'del': 0, 'mates_disagree': 0}
    weighted = dict([(k, 0) for k in counts])
    bad_reads = []
    unique_reads = set([])
    valid_reads = []
    
    reads = {}
    # get all mates
    count_ = 0
    for read in col.pileups:
        if read.alignment.flag & 0x4:
            continue
        unique_reads.add(read.alignment.qname)        
        reads[read.alignment.qname] = reads.get(read.alignment.qname, []) + [read]
    
    # check positions        
    for qname, pair in reads.items():
        base1, base2 = pair[0].alignment.seq[pair[0].qpos], None
        qual1, qual2 = ord(pair[0].alignment.qual[pair[0].qpos]) - 33, None
        isdel1, isdel2 = pair[0].is_del == 1, None        
        
        if len(pair) == 1:
            base, qual, isdel = base1, qual1, isdel1            
        else:
            base2 = pair[1].alignment.seq[pair[1].qpos]
            qual2 = ord(pair[1].alignment.qual[pair[1].qpos]) - 33
            isdel2 = pair[1].is_del == 1
            
            if isdel1 and isdel2:
                bad_reads.append(('del', qname))
                weighted['del'] += 1.0 
                counts['del'] += 1
                continue
            elif isdel1:
                base, qual, isdel = base2, qual2, False
            elif isdel2:
                base, qual, isdel = base1, qual1, False
            else:
                if base1 == base2:
                    # mates agree => check base  
                    base, qual, isdel = base1, max(qual1, qual2), False
                else:
                    # the mates do not agree at this position and have equal quality => discard read
                    weighted['mates_disagree'] += 1.0
                    counts['mates_disagree'] += 1
                    bad_reads.append(('mates_disagree', qname))
                    continue # <- important!
                pass
        if isdel:
            bad_reads.append(('del', qname))
            weighted['del'] += 1.0
            counts['del'] += 1            
        elif qual > cutoff and base != 'N':
            # print 'x', base, qual
            weighted[base] += 1.0
            counts[base] += 1
            valid_reads.append((qname, base))
        elif qual < cutoff:
            bad_reads.append(('low_qual', qname))
            weighted['low_qual'] += 1.0
            counts['low_qual'] += 1
        else:
            bad_reads.append(('N', qname))
            weighted['N'] += 1.0  
            counts['N'] += 1
        pass
    
    counts['T'] += counts['U']
    del counts['U']
    weighted['T'] += weighted['U']
    del weighted['U']
    counts['unique_reads'] = len(unique_reads)

    return counts, weighted, valid_reads

def count_reads_per_gene(rpg_raw, logfile=sys.stdout):
    # rpg_raw is dict with key=AGI, value=reads and reads is a list of tuples (id, AGI, bool:read_supports_mutation)
    # this function removes reads with conflicting SNPs from the count
    # in addition to checking locally (i.e. within AGI-genes), it checks reads mapping to different genes
    # in order to cover possible region-overlaps (which should occur only rarely)   
    rpg_filtered = {}
    all_reads = {}
    # flag all reads whether they are found supporting col (1), ped(2), or cause a conflict
    # because they support both col and ped at different positions (1 | 2 = 3) 
    for agi, reads in rpg_raw.items():
        for read in reads:
            # False = 0 + 1 => 1 (col-flag), True = 1 + 1 => 2 (ped-flag)            
            all_reads[read[0]] = all_reads.get(read[0], 0) | (int(read[-1]) + 1)
    # now remove all occurrences of conflicting reads
    # as well as reads that map different SNPs within the same gene (by making sets) 
    for agi, reads in rpg_raw.items():
        
        rpg_filtered[agi] = {'col': len(set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 1])),
                             'ped': len(set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 2])),
                             'conflict': len(set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 3]))}
        """
        Just for testing
        rpg_filtered[agi] = {'col': set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 1]),
                             'ped': set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 2]),
                             'conflict': set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 3])}
        """
        # rpg_filtered[agi]['checksum'] = len(rpg_raw[agi]) - sum(rpg_filtered[agi].values())
    
    #for agi, reads in rpg_raw.items():
    #    # for each gene/region of interest
    #    # 0. exclude reads that disagree in supporting mutation/wildtype at different snp sites
    #    support_ped = set([read[:-1] for read in reads if read[-1]])
    #    support_col = set([read[:-1] for read in reads if not read[-1]])
    #    undecided = support_col.intersection(support_ped)
    #            
    #    #rpg_filtered[agi] = {'ped': support_ped.difference(undecided), 
    #   #                     'col': support_col.difference(undecided),
    #    #                     'undecided': undecided}
    return rpg_filtered    
            
    
        
def analyse_snps(bamfile, snpfile_open, hit_mode, mult_counts=None, out=sys.stdout):
    header = ['', 'Pos_SNP', '#UniqueReads', '#Hits', '#Reads_Col', '#Reads_Ped', '#Score_Col', '#Score_Ped',
              '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']
    if hit_mode == 'intra':
        header[0] = 'AGI'
    elif hit_mode == 'inter':
        header[0] = 'Contig'
    else:
        sys.stderr.write('Wrong mode: %s. Exiting.\n' % hit_mode)
        sys.exit(1)
    out.write(','.join(map(str, header)) + '\n')
    
    analysed_snps = []    
    reads_per_gene = {}
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
                base_count, weighted, valid_reads = count_bases(col, mult_counts=mult_counts) 
    
        if base_count is not None: 
            
            
                       
            row = [None, end, 
                   base_count['unique_reads'],
                   sum(base_count.values()) - (base_count['unique_reads'] + base_count['del']),                   
                   base_count[attributes['refbase']],
                   base_count[attributes['mutation']],
                   weighted[attributes['refbase']],
                   weighted[attributes['mutation']],
                   base_count['N'],
                   base_count['low_qual'],
                   base_count['del'],
                   sum([base_count[v] for v in 'ACGT']) - \
                   (base_count[attributes['refbase']] + base_count[attributes['mutation']])]
            if hit_mode == 'intra':
                row[0] = attributes['gene_ID']
                rpg_key = row[0]                                                              
            elif hit_mode == 'inter':
                row[0] = contig
                rpg_key = '%s:%i' % (contig, end)                
            else:
                sys.stderr.write('Wrong mode: %s. Exiting.\n' % hit_mode)
                sys.exit(1)
            
            out.write(','.join(map(str, row)) + '\n')   
            analysed_snps.append(row)
            
            # list-set-list asserts that reads that cover multiple snps within the same gene
            # are only counted once (unless they disagree at different snp-sites, 
            # which is a different problem) 
            valid_reads = list(set([(read, rpg_key, base==attributes['mutation'])
                                    for read, base in valid_reads            
                                    if base in (attributes['refbase'], attributes['mutation'])]))
            if len(valid_reads) == 0:
                continue
                
            if rpg_key in reads_per_gene:                
                # reads_per_gene[attributes['gene_ID']] = reads_per_gene[attributes['gene_ID']].union(valid_reads)
                reads_per_gene[rpg_key].extend(valid_reads)
            else:
                reads_per_gene[rpg_key] = valid_reads
            
                         
            pass
        pass
    
    reads_per_gene = count_reads_per_gene(reads_per_gene)
    
    # if hit_mode == 'intra':         
    pickle.dump(reads_per_gene, open(bamfile.filename + '.' + hit_mode + '.rpg.dat', 'wb'))
    bamfile.close()
    
    # return analysed_snps, reads_per_gene
    pass
    


def main(argv):
    
    snp_fn = argv[0]
    if 'intra' in snp_fn:
        hit_mode = 'intra'
    elif 'inter' in snp_fn:
        hit_mode = 'inter'
    else:
        sys.stderr.write('Cannot find hit mode: %s. Exiting.\n' % snp_fn)
        sys.exit(1)
    
    # multi_hits = 'mult' in snp_fn
    multi_hits = False
    
    # bam_fn = snp_fn[:snp_fn.find('_')] + '.bam'
    bam_fn = argv[1]
    mult_counts = None
    if multi_hits:
        items = [line.strip().split('\t') 
                 for line in open(bam_fn + '.mult_counts', 'rb').readlines()]
        mult_counts = dict([(item[0], float(item[1])) for item in items])
        pass
    
    analyse_snps(pysam.Samfile(bam_fn, 'rb'), open(snp_fn, 'rb'), hit_mode, mult_counts=mult_counts, out=sys.stdout)
    pass

if __name__ == '__main__': main(sys.argv[1:])

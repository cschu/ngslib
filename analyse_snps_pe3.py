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

GENOME_PATH='/home/schudoma/ngslib/genomes/a_thaliana_TAIR10'

from count_bases import count_bases_old as count_bases

def count_reads_per_gene(rpg_raw, logfile=sys.stdout):
    """
    rpg_raw is dict with key=AGI, value=reads and reads is a list of tuples (id, AGI, bool:read_supports_mutation)
    this function removes reads with conflicting SNPs from the count
    in addition to checking locally (i.e. within AGI-genes), it checks reads mapping to different genes
    in order to cover possible region-overlaps (which should occur only rarely)
    """   
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
        
        col_set = [read for read in rpg_raw[agi] if all_reads[read[0]] == 1]
        ped_set = [read for read in rpg_raw[agi] if all_reads[read[0]] == 2]
        conflicts = [read for read in rpg_raw[agi] if all_reads[read[0]] == 3]
        
        
        
        rpg_filtered[agi] = {'col': len(set([read[0] for read in col_set])), 
                             'ped': len(set([read[0] for read in ped_set])), 
                             'conflict': len(set([read[0] for read in conflicts]))}
        # Just for testing
        # rpg_filtered[agi] = {'col': set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 1]),
        #                     'ped': set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 2]),
        #                     'conflict': set([read[0] for read in rpg_raw[agi] if all_reads[read[0]] == 3])}        
    return rpg_filtered    
            
def read_fasta(open_fn):
    return (''.join([line.strip() for line in open_fn.readlines()[1:]]))

def compute_consensus(seqs):
    consensus = ''
    if len(seqs) == 0:
        return consensus 
    for i, c in enumerate(seqs[0]):
        col = ''.join([seq[i] for seq in seqs if seq[i] != '-'])
        if len(col) == 0:
            consensus += '-'
        else:
            counts = sorted([(col.count(c), c) for c in set(list(col))], reverse=True)
            consensus += counts[0][1]
    return consensus


def get_primer_sequence(seqdata, offset=10):
    col, ped = [], []
    for seq, pos, is_mutant in seqdata:
        head = seq[max(0, pos-offset):pos]
        tail = seq[pos + 1:pos + 1 + offset]        
        primer = '-' * (offset - len(head)) + head + seq[pos] + tail + '-' * (offset - len(tail))        
        if is_mutant:
            ped.append(primer)
        else:
            col.append(primer)
    return compute_consensus(col), compute_consensus(ped)
        


def get_best_read_sequence(seqdata, offset=10):
    col = []
    ped = []
    for seq, pos, is_mutant in seqdata:
        d = (pos, len(seq) - pos)
        
        seq = seq[max(0, pos-offset):pos] + '*' + seq[pos] + '*' + seq[pos + 1:pos + 1 + offset]
        
        dat = (d, seq)
        if is_mutant:
            ped.append(dat)
        else:
            col.append(dat)
    col = [(None, None)] + sorted(col)
    ped = [(None, None)] + sorted(ped)
    
    # return col[-1], ped[-1]
    return col[-1][1], ped[-1][1]    
        
        
def analyse_snps(bamfile, snpfile_open, hit_mode, refpath=None, mult_counts=None, with_seqs=False, out=sys.stdout):
    header = ['', 'Pos_SNP', '#UniqueReads', '#Hits', '#Reads_Col', '#Reads_Ped', '#Score_Col', '#Score_Ped',
              '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']
    if hit_mode == 'intra':
        header[0] = 'AGI'
    elif hit_mode == 'inter':
        header[0] = 'Contig'
    else:
        sys.stderr.write('Wrong mode: %s. Exiting.\n' % hit_mode)
        sys.exit(1)
        
    if refpath is not None:        
        header.extend(['#match', '#bases', 'MAP_ERR'])
    
    if with_seqs:
        header.extend(['primer_col', 'primer_ped'])
    
    out.write(','.join(map(str, header)) + '\n')
        
    reads_per_gene = {}    
    contig = None
    for line in snpfile_open:
        snp = line.strip().split('\t')
        attributes = dict([attr.split('=') for attr in snp[8].split(';')])
        # next two if-conditions should always be False with new snp data
        # schedule for removal?
        if attributes['refbase'] not in 'ACGTU':
            continue
        if attributes['mutation'] not in 'ACGTU':
            continue
        
        refseq = None
        if refpath is not None:
            if contig is None or snp[0] != contig:            
                ref_fn = os.path.join(refpath, 'TAIR10_chr%c.fas' % snp[0][-1])
                refseq = read_fasta(open(ref_fn, 'rb'))
            
        
        contig, start, end = snp[0], int(snp[3]), int(snp[4])
        region = '%s:%i-%i' % (contig, start, end) 
        
        base_count = None
        flank_counts = []
        for col in bamfile.pileup(region=region):
            if start == col.pos:
                base_count, weighted, valid_reads = count_bases(col, mult_counts=mult_counts)
            elif refseq is not None:
                count, dummy, valid_reads_flank = count_bases(col, mult_counts=mult_counts)                
                try:
                    # this needs to be in a try-block because there might be ambiguity codes... ><;
                    flank_counts.append((count[refseq[col.pos]], len(valid_reads_flank)))
                except:
                    pass
                pass
    
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
            
            if refpath is not None:
                n_matches, n_bases = reduce(lambda x,y: (x[0]+y[0], x[1]+y[1]), flank_counts)
                mapping_error = 1.0 - (float(n_matches) / n_bases)
                row.extend([n_matches, n_bases, mapping_error])
            
            if hit_mode == 'intra':
                row[0] = attributes['gene_ID']
                rpg_key = row[0]                                                              
            elif hit_mode == 'inter':
                row[0] = contig
                rpg_key = '%s:%i' % (contig, end)                
            else:
                sys.stderr.write('Wrong mode: %s. Exiting.\n' % hit_mode)
                sys.exit(1)
            
            # list-set-list asserts that reads that cover multiple snps within the same gene
            # are only counted once (unless they disagree at different snp-sites, 
            # which is a different problem) 
            
            seqs = [(seq, qpos, base==attributes['mutation']) 
                    for read, base, seq, qpos in valid_reads
                    if base in (attributes['refbase'], attributes['mutation'])
                    and seq is not None]
            # colseq, pedseq = get_best_read_sequence(seqs)
            if with_seqs:
                colseq, pedseq = get_primer_sequence(seqs)
                row.extend([colseq, pedseq])
            out.write(','.join(map(str, row)) + '\n')
            
            valid_reads = list(set([(read, rpg_key, base==attributes['mutation'])
                                    for read, base, seq, qpos in valid_reads            
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
    
    pickle.dump(reads_per_gene, open(os.path.basename(bamfile.filename) + '.' + hit_mode + '.rpg.full.dat', 'wb'))
    
    reads_per_gene = count_reads_per_gene(reads_per_gene)
    
    # if hit_mode == 'intra':         
    pickle.dump(reads_per_gene, open(os.path.basename(bamfile.filename) + '.' + hit_mode + '.rpg.dat', 'wb'))
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
    
    bam_fn = argv[1]
    
    fo = os.path.basename(bam_fn)
    version = '003'
    fo = version + '.' + fo + '.' + hit_mode + '.' + 'snps.csv'
    refpath = None
    analyse_snps(pysam.Samfile(bam_fn, 'rb'), open(snp_fn, 'rb'), hit_mode, refpath=refpath, out=open(fo, 'wb'))
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/python

import os
import sys
import math

import pysam

from snptools import SNPDict

# Peol-shoot_TGACCA_L007_R1_001.paired.fq.gsnap.concordant_uniq.bam.samsnp

MIN_NREADS = 10
COL_HEADERS = ['Contig', 'Start', 'End', 'Gene', 'Genetype', 
               '#reads', '#snps', '#covered',
               '#support_ref', '#support_mut', 
               'undecided', 'fr_ref', 'fr_snp', 
               '#rsupport_ref', '#rsupport_mut', '#valid_reads', '#bad',
               'rfr_ref', 'rfr_mut', 'call',
               'comment']

class DataLine(dict):
    def __init__(self, keys=[]):
        for key in keys:
            self[key] = ''
        pass
    def set_range(self, keys=[], values=[]):
        indices = range(len(self), len(self) + len(keys) + 1)
        for k, v in zip(keys, zip(indices, map(str, values))):
            self[k] = v
        pass
    def set(self, key, value):
        self[key] = (len(self), str(value))
        pass
    def get_string(self, sep=','):
        return sep.join([x[1] for x in sorted(self.values(), key=lambda x:x[0])])
    def get_headerstring(self, sep=','):
        return sep.join([x[0] for x in sorted(self.items(), key=lambda x:x[1][0])])
    pass
    

#
def get_reads(samfile, contig, start, end):
    return [read.qname 
            for read in samfile.fetch(contig, start, end)]

#
def get_geneid(gffline):
    comments_d = parse_gff_comments(gffline[8])     
    return gffline[:1] + gffline[3:5] + [comments_d.get('ID', 'N/A'),
                                         comments_d.get('Note', 'N/A')]
#
def parse_gff_comments(string):
    return dict([x.split('=')
                 for x in string.strip().split(';')])

#
def read_polymorphs(open_gff):
    polymorphs = {}
    for gffline in open_gff:
        gffline = gffline.strip().split('\t')
        comments = parse_gff_comments(gffline[8])
        snp = (gffline[0], int(gffline[3]) - 1, 
               int(gffline[4]) - 1, gffline[6])
        gene_id = comments['ID']
        polymorphs[gene_id] = polymorphs.get(gene_id, []) + [snp]
    return polymorphs

#    
def call_snp(fr_ref, fr_mut, cutoff=0.75, guards=(0.25, 0.75)):
    if fr_ref == 'NaN' or fr_mut == 'NaN':
        return 'n/a'
    elif fr_ref + fr_mut > cutoff:
        if fr_ref/(fr_ref + fr_mut) > guards[1]:
            return 'ref'
        elif fr_ref/(fr_ref + fr_mut) < guards[0]:
            return 'mut'
        else:
            return 'undecided'
    else:
        return 'noisy'         

#
def count_bases(samfile, contig, pos, cutoff=-5):
    counts = {'bad': 0}
    for col in samfile.pileup(contig, start=pos, end=pos+1):
        if col.pos == pos:
            for read in col.pileups:
                if read.is_del == 0:
                    base = read.alignment.seq[read.qpos]
                    qual = ord(read.alignment.qual[read.qpos]) - 33
                    if qual > cutoff and base != 'N':
                        counts[base] = counts.get(base, 0) + 1
                    else:
                        counts['bad'] = counts.get('bad', 0) + 1
                else:
                    counts['bad'] = counts.get('bad', 0) + 1
    return counts

#    
def classify_position_from_basecount(basecount, refbase, snpbase, index):
    comment = ''    
    if len(basecount) == 0:
        comment = 'not-supported-%i' % index
    elif len(basecount) == 1:
        if (refbase in basecount) or (snpbase in basecount):
            comment = 'unique-snp-%i' % index
        else:
            comment = 'weird-unique-snp-%i' % index
    elif len(basecount) == 2:
        if (refbase in basecount) or (snpbase in basecount):
            comment = 'undecided-snp-%i' % index
        else:
            comment = 'weird-multisnp-%i' % index
    else:
        comment = 'multisnp-%i-%i' % (len(basecount), index)
        
    return comment

#    
def process_snps(snps, snp_d, contig, samfile):
    snp_classes = []
    snp_fracs = []
    snps_supporting_mutation = 0
    snps_supporting_reference = 0 
    covered_snps = []  
    
    reads_supporting_mutation = 0
    reads_supporting_reference = 0
    total_reads = 0    
    
    for i, snp in enumerate(snps):
        basecount = count_bases(samfile, contig, snp[1])
        n_reads = sum(basecount.values())
        
        snpline = snp_d.get((contig, snp[1] + 1))
        refbase = float(basecount.get(snpline.refbase, 0.0))
        snpbase = float(basecount.get(snpline.mutation, 0.0))
        
        reads_supporting_mutation += snpbase
        reads_supporting_reference += refbase
        total_reads += sum(basecount.values()) - basecount['bad']

        if n_reads > 0:
            snp_fracs.append((refbase/n_reads, snpbase/n_reads))
        else:
            snp_fracs.append(None)
            
        cur_comment = classify_position_from_basecount(basecount,
                                                       snpline.refbase, 
                                                       snpline.mutation, i)
        snp_classes.append(cur_comment)
        
        if len(basecount) == 0:
            pass
        elif refbase > snpbase:
            snps_supporting_reference += 1.0
        elif snpbase > refbase:
            snps_supporting_mutation += 1.0
        else:
            snps_supporting_reference += 0.5
            snps_supporting_mutation += 0.5
        
        # if len(snp_reads) > 0:
        if n_reads > 0:
            # covered_snps.append((snp, len(snp_reads)))
            covered_snps.append((snp, n_reads))
            
    counts = (snps_supporting_reference, snps_supporting_mutation,
              reads_supporting_reference, reads_supporting_mutation, total_reads, basecount['bad'])
    
    return snp_classes, snp_fracs, covered_snps, counts

#
def process_gff(open_gff, polymorphs, snp_d, samfile, fo, fo2, min_reads=10):
    
    header = False
    
    for gffline in open_gff:
        gffline = gffline.split('\t')
        start, end = map(lambda x:x-1, map(int, gffline[3:5]))
        contig = gffline[0]
        
        reads = get_reads(samfile, contig, start, end)
        # n_reads = len(reads)
        
        idstr = get_geneid(gffline)
        snps = polymorphs.get(idstr[3], [])
        
        if len(snps) > 0:
            classes, fracs, covered, counts = process_snps(snps, snp_d,
                                                           contig, samfile)
            undecided = len(snps) - sum(counts[:1])
            n = len(fracs)
            mean_fr_ref = sum([x[0] for x in fracs if not x is None]) / n
            mean_fr_snp = sum([x[1] for x in fracs if not x is None]) / n
            
            fr_ref, fr_mut = 'NaN', 'NaN'
            if counts[4] != 0.0:
                fr_ref, fr_mut = counts[2]/float(counts[4]), counts[3]/float(counts[4])

            if counts[4] >= min_reads:
                data = DataLine()
                data.set_range(keys=COL_HEADERS[:5], values=idstr)
                data.set_range(keys=COL_HEADERS[5:8], values=[counts[4], len(snps), len(covered)])
                data.set_range(keys=['#support_ref', '#support_mut', '#valid_reads', '#bad'], values=counts[2:])
                data.set_range(keys=['fr_ref', 'fr_mut', 'call'], values=[fr_ref, fr_mut, call_snp(fr_ref, fr_mut)])
                data.set('#undecided', counts[4] - (fr_ref + fr_mut))
                data.set('comment', '|'.join(classes))
                
                """
                data = [counts[4], len(snps), len(covered), counts[0], counts[1], undecided,
                        mean_fr_ref, mean_fr_snp,
                        counts[2], counts[3], counts[4], counts[5],
                        fr_ref, fr_mut,
                        call_snp(fr_ref, fr_mut), 
                        '|'.join(classes)]
                """
                
                # outstr = '%s\n' % '\t'.join(idstr + map(str, data))
                outstr = data.get_string(sep='\t')
                if not header:
                    header = True
                    fo.write(data.get_headerstring(sep='\t'))
                fo.write('%s\n' % outstr)
                fo2.write('%s\n' % outstr)
                for read in reads:
                    fo2.write('%s\n' % read)
    
    return None

###
def main(argv):

    samfile = pysam.Samfile(argv[0], 'rb')
    fo = sys.stdout
    fo = open('%s.covered_genes.csv' % argv[0].rstrip('.bam'), 'w')    
    fo2 = open('%s.covered_genes_with_reads.csv' % argv[0].rstrip('.bam'), 'w')

    # fo.write('%s\n' % ','.join(COL_HEADERS))

    polymorphs = read_polymorphs(open(argv[1]))
    snp_d = SNPDict(open(argv[2]))    
    process_gff(open(argv[3]), polymorphs, snp_d, samfile, fo, fo2, min_reads=MIN_NREADS)
    fo2.close()
    fo.close()
    
    return None
    
if __name__ == '__main__': main(sys.argv[1:])

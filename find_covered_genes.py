#!/usr/bin/python

import os
import sys
import math

import pysam

from snptools import SNPDict


# Peol-shoot_TGACCA_L007_R1_001.paired.fq.gsnap.concordant_uniq.bam.samsnp

def parse_gff_comments(string):
    return dict([x.split('=')
                 for x in string.strip().split(';')])


###
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

###
def count_bases(samfile, contig, pos, cutoff=-5):
    # print contig, pos
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
    del counts['bad']
    return counts
    
    
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


def get_reads(samfile, contig, start, end):
    return [read.qname 
            for read in samfile.fetch(contig, start, end)]
    
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
        snp_reads = [read
                     for read in samfile.fetch(contig, 
                                               start=snp[1], end=snp[1] + 1)]
        basecount = count_bases(samfile, contig, snp[1])
        snpline = snp_d.get((contig, snp[1] + 1))
        refbase = float(basecount.get(snpline.refbase, 0.0))
        snpbase = float(basecount.get(snpline.mutation, 0.0))
        
        reads_supporting_mutation += snpbase
        reads_supporting_reference += refbase
        total_reads += sum(basecount.values())
        
        if len(basecount) > 0:
            snp_fracs.append((refbase/sum(basecount.values()),
                              snpbase/sum(basecount.values())))
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
        
        if len(snp_reads) > 0:
            covered_snps.append((snp, len(snp_reads)))
            
    counts = (snps_supporting_reference, snps_supporting_mutation,
              reads_supporting_reference, reads_supporting_mutation, total_reads)
    return snp_classes, snp_fracs, covered_snps, counts
        
    

def process_gff(open_gff, polymorphs, snp_d, samfile, fo, fo2, min_reads=10):
    
    for gffline in open_gff:
        gffline = gffline.split('\t')
        start, end = map(lambda x:x-1, map(int, gffline[3:5]))
        contig = gffline[0]
        reads = get_reads(samfile, contig, start, end)
        n_reads = len(reads)
        
        if n_reads >= min_reads:
            comments_d = parse_gff_comments(gffline[8])
            geneid = comments_d.get('ID', 'N/A')
            idstr = gffline[:1] + gffline[3:5] + [geneid, 
                                                  comments_d.get('Note', 'N/A')]
            snps = polymorphs.get(geneid, [])
            if len(snps) == 0:
                continue
            
            classes, fracs, covered, counts = process_snps(snps, snp_d, 
                                                           contig, samfile)
            undecided = len(snps) - sum(counts)
            n = len(fracs)
            mean_fr_ref = sum([x[0] for x in fracs if not x is None]) / n
            mean_fr_snp = sum([x[1] for x in fracs if not x is None]) / n
            
            fr_ref, fr_mut = 'NaN', 'NaN'
            if counts[4] != 0.0:
                fr_ref, fr_mut = counts[2]/float(counts[4]), counts[3]/float(counts[4])
                
            
            data = [n_reads, len(snps), len(covered), counts[0], counts[1], undecided,
                    mean_fr_ref, mean_fr_snp,
                    counts[2], counts[3], counts[4],
                    fr_ref, fr_mut, 
                    '|'.join(classes)]
            
            outstr = '%s\n' % '\t'.join(idstr + map(str, data))
            fo.write(outstr)
            fo2.write(outstr)
            for read in reads:
                fo2.write('%s\n' % read)
    
    return None
    
        
        


###
def main(argv):

    MIN_NREADS = 10
    col_headers = ['Contig', 'Start', 'End', 'Gene', 'Genetype', 
                   '#reads', '#snps', '#covered',
                   '#support_ref', '#support_mut', 
                   'undecided', 'fr_ref', 'fr_snp', 
                   '#rsupport_ref', '#rsupport_mut', '#valid_reads',
                   'rfr_ref', 'rfr_mut', 
                   'comment']


    samfile = pysam.Samfile(argv[0], 'rb')
    fo = open('%s.covered_genes.csv' % argv[0].rstrip('.bam'), 'w')
    fo2 = open('%s.covered_genes_with_reads.csv' % argv[0].rstrip('.bam'), 'w')

    fo.write('%s\n' % ','.join(col_headers))

    polymorphs = read_polymorphs(open(argv[1]))
    snp_d = SNPDict(open(argv[2]))    
    
    process_gff(open(argv[3]), polymorphs, snp_d, samfile, fo, fo2, min_reads=MIN_NREADS)
    fo2.close()
    fo.close()
    
    
    
    return None
    for gffline in open(argv[3]):
        gffline = gffline.split('\t')
        # print gffline
        start, end = map(lambda x:x-1, map(int, gffline[3:5]))
        contig = gffline[0]
        # print contig, start, end

        reads = []        
        for read in samfile.fetch(contig, start, end):
            # print read.pos,
            reads.append(read.qname)
        # print 
        if len(reads) >= MIN_NREADS:    
            comments_d = parse_gff_comments(gffline[8])
            geneid = comments_d.get('ID', 'N/A')

            idstr = [geneid,
                     comments_d.get('Note', 'N/A')]
            idstr = gffline[0:1] + gffline[3:5] + idstr        

            snps = polymorphs.get(geneid, [])
            if len(snps) == 0:
                continue
            
            comment = []
            covered_snps = []
            snps_supporting_mutation = 0
            snps_supporting_reference = 0
            
            snp_fracs = []
            for i, snp in enumerate(snps):
                # print snp                
                snp_reads = [read 
                             for read in samfile.fetch(contig, 
                                                       start=snp[1],
                                                       end=snp[1]+1)]
                basecount = count_bases(samfile, contig, snp[1])                
                # print basecount
                snpline = snp_d.get((contig, snp[1] + 1))
                # print snpline
                # maxbase = max(basecount.values())
                refbase = float(basecount.get(snpline.refbase, 0.0))
                snpbase = float(basecount.get(snpline.mutation, 0.0))
                if len(basecount) > 0:
                    snp_fracs.append((refbase/sum(basecount.values()),
                                      snpbase/sum(basecount.values())))
                else:
                    snp_fracs.append(None)
                
                cur_comment = classify_position_from_basecount(basecount, 
                                                               snpline.refbase, 
                                                               snpline.mutation, i)
                comment.append(cur_comment)

                if len(basecount) == 0:
                    pass
                elif refbase > snpbase:
                    snps_supporting_reference += 1
                elif snpbase > refbase:
                    snps_supporting_mutation += 1
                else:
                    snps_supporting_reference += 0.5
                    snps_supporting_mutation += 0.5
                
                # sys.exit(1)

                # DEBUG BLOCK
                #snp_reads = [col for col in samfile.pileup(contig,
                #                                           snp[1])]
                # print sorted([int(x.pos) for x in snp_reads])[:30]
                # print sorted([int(x.pos) for x in snp_reads])[-30:]
                if len(snp_reads) > 0:
                    covered_snps.append((snp, len(snp_reads)))
            
            undecided = len(snps)
            undecided -= (snps_supporting_reference + snps_supporting_mutation)
            # print covered_snps, len(covered_snps)  
            
            n = len(snp_fracs)
            mean_fr_ref = sum([x[0] 
                               for x in snp_fracs if not x is None]) / n
            mean_fr_snp = sum([x[1] 
                               for x in snp_fracs if not x is None]) / n
            
  
            """
            try:
                max_reads_covering_snp = max([x[1] for x in covered_snps])
            except:
                max_reads_covering_snp = 0
            """
            data = [len(reads), len(snps), len(covered_snps),
                    snps_supporting_reference,
                    snps_supporting_mutation,
                    undecided,
                    mean_fr_ref,
                    mean_fr_snp,
                    '|'.join(comment)]

                    
            outstr = '%s\n' % '\t'.join(idstr + map(str, data))
                       
            fo.write(outstr)
            fo2.write(outstr)
            for read in reads:
                fo2.write('%s\n' % read)

            # ID=AT1G01040;Note=protein_coding_gene;Name=AT1G01040
        pass
    fo2.close()
    fo.close()

    return None

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/python

import os
import sys
import math

import pysam

from snptools import SNPDict
import gff_helpers

# Peol-shoot_TGACCA_L007_R1_001.paired.fq.gsnap.concordant_uniq.bam.samsnp

MIN_NREADS = 3
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
    return [(read.qname, read.tags) 
            for read in samfile.fetch(contig, start, end)]


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
    
INTERVAL_GUARDS = map(lambda x:x/100., range(10,105,10))    
    
def call_snp_interval(fr_ref, guards=INTERVAL_GUARDS):
    for guard in guards:
        if guard >= fr_ref:
            return (guards.index(guard) + 1)/float(len(guards))
    return -1
#
def check_nomapping(read):
    return 'nomapping' in dict(read.tags).get('RG', '') 

def count_bases(samfile, contig, pos, cutoff=-5):
    counts = {'bad': 0, 'nomap': 0}
    # print 'COUNT_BASES:', contig, pos, pos+1
    
    """for col in samfile.pileup(contig, start=pos, end=pos+1):
        print col
    """
    for col in samfile.pileup(contig, start=pos, end=pos+1, stepper='all'):
        
        if col.pos == pos:
            readcount = 0                                    
            for read in col.pileups:
                readcount += 1
                #if check_nomapping(read):
                #    counts['nomap'] = counts.get('nomap', 0) + 1
                #elif read.is_del == 0:
                if read.is_del == 0:
                    base = read.alignment.seq[read.qpos]
                    qual = ord(read.alignment.qual[read.qpos]) - 33
                    if qual > cutoff and base != 'N':
                        counts[base] = counts.get(base, 0) + 1
                    else:
                        counts['bad'] = counts.get('bad', 0) + 1
                else:
                    counts['bad'] = counts.get('bad', 0) + 1
            # print 'READCOUNT:', readcount
    return counts

#    
def classify_position_from_basecount(basecount, refbase, snpbase, index):
    bad_count = basecount['bad']
    nomap_count = basecount['nomap']
    del basecount['bad']
    del basecount['nomap']
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
    # print basecount, refbase, snpbase, comment
    basecount['bad'] = bad_count    
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

    bad_count = 0

    actual_snps = len(snps)    
    for i, snp in enumerate(snps):
        snpline = snp_d.get((contig, snp[1] + 1))
        # print 'SNP', (contig, snp[1] + 1), snpline, snp_d.get((contig, snp[1])), snp_d.get((contig, snp[1] - 1))
        if snpline is None:
            actual_snps -= 1
            continue
        
        basecount = count_bases(samfile, contig, snp[1])
        # print basecount
        # sys.exit(1)
        bad_count += basecount['bad']
        n_reads = sum(basecount.values())
        
         
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
              reads_supporting_reference, reads_supporting_mutation, total_reads, bad_count)
    
    return snp_classes, snp_fracs, covered_snps, counts, actual_snps

#
def process_gff(open_gff, polymorphs, snp_d, samfile, fo, fo2, min_reads=10):
    
    header = False
    
    total_gfflines = 0
    has_snps = 0
    enough_reads = 0 
    
    for gffline in open_gff:
        total_gfflines += 1
        gffline = gffline.split('\t')
        start, end = map(lambda x:x-1, map(int, gffline[3:5]))
        contig = gffline[0]
        
        # reads = get_reads(samfile, contig, start, end)
        # print reads[0]
        # sys.exit(1)
        # n_reads = len(reads)
        
        idstr = gff_helpers.get_geneid(gffline)
        snps = polymorphs.get(idstr[3], [])
        # print idstr#, snps
         
        if len(snps) > 0:
            has_snps += 1
            # actual snps variable is needed due to polymorph-library containing 'indel-snps' 
            classes, fracs, covered, counts, actual_snps = process_snps(snps, snp_d,
                                                                        contig, samfile)            
            
            fr_ref, fr_mut = 'NaN', 'NaN'
            if counts[4] != 0.0:
                fr_ref, fr_mut = counts[2]/float(counts[4]), counts[3]/float(counts[4])

            count4_gt_min_reads = False
            if counts[4] >= min_reads:
                enough_reads += 1
                count4_gt_min_reads = True
            # if True:
                data = DataLine()
                data.set_range(keys=COL_HEADERS[:5], values=idstr)
                data.set_range(keys=COL_HEADERS[5:8], values=[counts[4], actual_snps, len(covered)])
                data.set_range(keys=['#support_ref', '#support_mut'], values=counts[2:4])
                data.set('#undecided', counts[4] - (counts[2] + counts[3]))
                data.set_range(keys=['#total_reads', '#bad'], values=counts[4:])
                data.set_range(keys=['fr_ref', 'fr_mut', 'call'], values=[fr_ref, fr_mut, call_snp_interval(fr_ref)])#, fr_mut)])
                
                data.set('count4', str(count4_gt_min_reads))
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
                    fo.write('%s\n' % data.get_headerstring(sep='\t'))
                fo.write('%s\n' % outstr)
                #fo2.write('%s\n' % outstr)
                #for read in reads:
                #    fo2.write('%s\n' % read[0])
    print total_gfflines, has_snps, enough_reads 
    
    return None

###
def main(argv):

    samfile = pysam.Samfile(argv[0], 'rb')
    fo = sys.stdout
    fo = open('%s.covered_genes.csv' % argv[0].rstrip('.bam'), 'w')    
    # fo2 = open('%s.covered_genes_with_reads.csv' % argv[0].rstrip('.bam'), 'w')
    fo2=sys.stdout

    # fo.write('%s\n' % ','.join(COL_HEADERS))
    
    # tair10_ped-0_polymorphic_genes.gff
    polymorphs = gff_helpers.read_polymorphs(open(argv[1]))
    # ped-0-snps_no-indels.txt
    snp_d = SNPDict(open(argv[2]))  
    # tair10_genes.gff
    process_gff(open(argv[3]), polymorphs, snp_d, samfile, fo, fo2, min_reads=MIN_NREADS)
    # fo2.close()
    fo.close()
    samfile.close()
    return None
    
if __name__ == '__main__': main(sys.argv[1:])

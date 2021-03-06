#!/usr/bin/env python
'''
Created on Nov 16, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re

import pysam


def get_increment(mult_counts, qname):
    if '_' in qname:
        qname = qname[:qname.rfind('_')]
    try:
        return 1.0 / mult_counts[qname]
    except:
        return 1.0
    pass

def find_unique_mate_pairs(reads):
    
    first = [read for read in reads 
             if (read.alignment.flag & 0x40 or read.alignment.flag & 0x1 == 0)]
    second = [read for read in reads if read.alignment.flag & 0x80]
    # print first
    # print second

    unique_pairs = []
    for k, read1 in enumerate(first):
        ra1 = read1.alignment
        mate, j = None, -1
        for i, read2 in enumerate(second):
            ra2 = read2.alignment
            if ra1.pnext == ra2.pos and ra1.pos == ra2.pnext:
                mate, j = read2, i
                break
        new_pair = [read1]
        if mate is not None:
            new_pair.append(mate)            
            del second[j]
        unique_pairs.append(new_pair)
    for read2 in second:
        unique_pairs.append([read2])
    
    qname = reads[0].alignment.qname
    return zip([qname + '_%i' % i for i in xrange(len(unique_pairs))], unique_pairs)



def count_bases(col, cutoff=-5, mult_counts=None):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 
              'low_qual': 0, 'del': 0, 'mates_disagree': 0}
    weighted = dict([(k, 0) for k in counts])
    bad_reads = []
    unique_reads = set([])
    
    reads = {}
    # get all mates and multimappings
    count_ = 0
    for read in col.pileups:
        if read.alignment.flag & 0x4:
            continue
        unique_reads.add(read.alignment.qname)        
        reads[read.alignment.qname] = reads.get(read.alignment.qname, []) + [read]
    # print unique_reads, len(unique_reads)
    
    # for read in col.pileups:
    #    count_ += 1
    #    if read.alignment.flag & 0x40:
    #        print 'xx'
    #    print 'dINC', get_increment(mult_counts, read.alignment.qname)
    #    unique_reads.add(read.alignment.qname)        
    #    if read.is_del == 0:
    #        reads[read.alignment.qname] = reads.get(read.alignment.qname, []) + [read]
    #        
    #    else:
    #        bad_reads.append(('del', read.alignment.qname))
    #        weighted['del'] += get_increment(mult_counts, read.alignment.qname)
    #        counts['del'] += 1
    #print reads
    #print count_
        
    # sort out the multimappings in paired-end reads
    if mult_counts is not None:
        tmp_reads = []    
        for qname in reads:
            tmp_reads += find_unique_mate_pairs(reads[qname]) 
            pass
        reads = dict(tmp_reads)
    # print reads
    
    # check positions        
    for qname, pair in reads.items():
        # print 'INC', get_increment(mult_counts, qname)
        # print qname, pair, len(pair)
        print qname

        base1, base2 = pair[0].alignment.seq[pair[0].qpos], None
        qual1, qual2 = ord(pair[0].alignment.qual[pair[0].qpos]) - 33, None
        isdel1, isdel2 = pair[0].is_del == 1, None 
        
        
        if len(pair) == 1:
            # print 'x', base1, qual1, isdel1
            base, qual, isdel = base1, qual1, isdel1            
        else:
            base2 = pair[1].alignment.seq[pair[1].qpos]
            qual2 = ord(pair[1].alignment.qual[pair[1].qpos]) - 33
            isdel2 = pair[1].is_del == 1
            
            if isdel1 and isdel2:
                bad_reads.append(('del', qname))
                weighted['del'] += get_increment(mult_counts, qname)
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
                # elif qual1 == qual2:
                else:
                    # the mates do not agree at this position and have equal quality => discard read
                    weighted['mates_disagree'] += get_increment(mult_counts, qname) 
                    counts['mates_disagree'] += 1
                    bad_reads.append(('mates_disagree', qname))
                    continue # <- important!
                # elif qual1 > qual2:
                #    # mates disagree with mate1 having higher quality
                #    base, qual = base1, qual1
                # else:
                #    # mates disagree with mate2 having higher quality
                #    base, qual = base2, qual2
            
        if isdel:
            bad_reads.append(('del', qname))
            weighted['del'] += get_increment(mult_counts, qname)
            counts['del'] += 1            
        elif qual > cutoff and base != 'N':
            print 'x', base, qual
            weighted[base] += get_increment(mult_counts, qname)
            counts[base] += 1
        elif qual < cutoff:
            bad_reads.append(('low_qual', qname))
            weighted['low_qual'] += get_increment(mult_counts, qname)
            counts['low_qual'] += 1
        else:
            bad_reads.append(('N', qname))
            weighted['N'] += get_increment(mult_counts, qname)  
            counts['N'] += 1
        # print counts
        pass
    
    counts['T'] += counts['U']
    del counts['U']
    weighted['T'] += weighted['U']
    del weighted['U']
    
    counts['unique_reads'] = len(unique_reads)
    print counts
    # print counts
    # print weighted

    return counts, weighted
        
        
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
                base_count, weighted = count_bases(col, mult_counts=mult_counts) 
    
        if base_count is not None:            
            out = [None, end, 
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
    if 'intra' in snp_fn:
        hit_mode = 'intra'
    elif 'inter' in snp_fn:
        hit_mode = 'inter'
    else:
        sys.stderr.write('Cannot find hit mode: %s. Exiting.\n' % snp_fn)
        sys.exit(1)
    
    multi_hits = 'mult' in snp_fn
    
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

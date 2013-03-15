#!/usr/bin/env python
'''
Created on Mar 14, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import pickle

import pysam

TRANSCRIPT_ATTR = ['contig', None, 'type', 'start', 'end', None, 'strand', None, 'id_']
TRANSCRIPT_CASTS = [str, None, str, int, int, None, str, None, str]

PILEUP_ATTR = ['support_refbase', 'support_mutation', 'N_at_pos', 'mismatches', 'lowqual', 'indel', 'n_aligned_reads']
R_PILEUP_ATTER = ['r_support_refbase', 'r_support_mutation', 'r_N_at_pos', 'r_mismatches', 'r_lowqual', 'r_indel', 'r_n_aligned_reads']

from count_bases import count_bases

class Transcript(object):
    def __init__(self, attributes):
        self.snps = []
        for name, val, cast in [x 
                                for x in zip(TRANSCRIPT_ATTR, attributes, TRANSCRIPT_CASTS) 
                                if not x[0] is None]:
            if name == 'id_':
                val = val.split(';')[0].lstrip('ID=')
                pass
            setattr(self, name, cast(val))
        pass
    def __str__(self):
        return '\t'.join(map(str, [self.id_, self.contig, self.type, self.start, self.end, self.strand]))
    pass
    
            
class SNP_Position(object):
    def __init__(self, attributes):        
        self.pos_0 = int(attributes[0])
        self.pos_1 = int(attributes[1])
        self.count = 0
        self.test = 0
        self.is_covered = False
        for name, val in [x.split('=') for x in attributes[2].split(';')]:
            try:
                val = float(val)
            except:
                pass
            setattr(self, name, val)
        pass
    def __str__(self):        
        if hasattr(self, 'support_refbase'):
            readstats = '\t'.join(map(str, [self.support_refbase, self.support_mutation,
                                            len(self.r_support_refbase), len(self.r_support_mutation)]))
        else:
            readstats = '0\t0\t0\t0'        
        return '\t'.join(map(str, [self.pos_0, self.pos_1, self.refbase, self.mutation, readstats]))
    def get_region(self, contig):
        return '%s:%i-%i' % (contig, self.pos_0, self.pos_1)
    def add_pileup(self, pileup_count, reads_at_pos):
        self.is_covered = True
        for attr in PILEUP_ATTR:
            if not hasattr(self, attr):
                setattr(self, attr, 0)
        for attr in R_PILEUP_ATTER:
            if not hasattr(self, attr):
                setattr(self, attr, [])
        
        self.support_refbase += pileup_count[self.refbase]
        self.support_mutation += pileup_count[self.mutation]
        self.N_at_pos += pileup_count['N']
        self.mismatches += sum([pileup_count[c] for c in 'ACGT']) - \
        (pileup_count[self.refbase] + pileup_count[self.mutation])
        self.lowqual += pileup_count['low_qual']
        self.indel += pileup_count['del']
        self.n_aligned_reads += pileup_count['unique_reads'] - self.indel
        
        self.r_support_refbase += reads_at_pos[self.refbase]
        self.r_support_mutation += reads_at_pos[self.mutation]
        self.r_N_at_pos += reads_at_pos['N']
        
        abc = list(set(list('ACGT')) - set([self.refbase, self.mutation]))
        for c in abc:
            self.r_mismatches += reads_at_pos[c]
        self.r_lowqual += reads_at_pos['low_qual']
        self.r_indel += reads_at_pos['del']
        self.r_n_aligned_reads += list(set(reads_at_pos['unique_reads']) - set(self.r_indel))        
        pass
    def cleanse_counts(self, valid_reads):        
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) - set(self.r_support_refbase))
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) - set(self.r_support_mutation))                
        self.r_support_refbase = list(set(self.r_support_refbase).intersection(valid_reads))
        self.r_support_mutation = list(set(self.r_support_mutation).intersection(valid_reads))        
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) | set(self.r_support_refbase))
        self.r_n_aligned_reads = list(set(self.r_n_aligned_reads) | set(self.r_support_mutation))
        pass
                
    pass
      

def read_transcript_data(open_fn):
    transcript_d = {}
    snp_d = {}
    for line in open_fn:
        line = line.strip().split('\t')
        transcript, snp = Transcript(line[:9]), SNP_Position(line[9:])
        
        snp_key = (transcript.contig, snp.pos_1)
        if snp_key in snp_d:
            snp = snp_d[snp_key][0]
        else:
            snp_d[snp_key] = [snp]
        snp_d[snp_key].append(transcript.id_)
        snp.count += 1    
                
        if transcript.id_ not in transcript_d:
            transcript_d[transcript.id_] = transcript
        transcript_d[transcript.id_].snps.append(snp)
        
    return transcript_d, snp_d
        
def filter_ambiguous_reads(reads):
    all_reads = {}
    # flag all reads whether they are found supporting col (1), ped(2), or cause a conflict
    # because they support both col and ped at different positions (1 | 2 = 3)
    for read in reads:         
        all_reads[read[0]] = all_reads.get(read[0], 0) | (int(read[-1]) + 1)
    # now remove all occurrences of conflicting reads
    # as well as reads that map different SNPs within the same gene (by making sets)
    return [read for read in all_reads if all_reads[read] != 3]
      
        
         
        
def process_pileups(bamfile, snp_d, read_checklist):
    N = len(snp_d)
    i = 1.0    
    for k in snp_d:
        if (i*100/N) % 5 == 0: 
            sys.stderr.write('%i SNPs processed (%i%%)\n' % (i, i/N * 100))
        i += 1.0
        snp = snp_d[k][0]
        snp.test += 1
        base_count, reads_at_pos = None, None        
        for col in bamfile.pileup(region=snp.get_region(k[0])):
            if snp.pos_0 == col.pos:
                base_count, weighted, valid_reads, reads_at_pos = count_bases(col)
                break
        if base_count is None or reads_at_pos is None:
            continue
        snp.add_pileup(base_count, reads_at_pos)
        
        valid_reads = set([(read, base==snp.mutation)
                           for read, base, seq, qpos in valid_reads            
                           if base in (snp.refbase, snp.mutation)])
        read_checklist.update(valid_reads)
    
    # This is a simulation to see how ambiguous SNPs can be removed from the read count
    # print snp_d.keys()
    # snp_d[('Chr4', 10952350)][0].r_support_mutation.append('HWI-ST377:127:D0PHGACXX:7:1306:11252:6941')
    # read_checklist.add(('HWI-ST377:127:D0PHGACXX:7:1306:11252:6941', True))
    # print 'x', sorted(snp_d[('Chr4', 10952350)][0].r_support_refbase[:10])
    # print 'x', sorted(snp_d[('Chr4', 10951240)][0].r_support_refbase[:10])
    
    filtered_reads = set(filter_ambiguous_reads(read_checklist))
    # print len(read_checklist), len(filtered_reads) 
    for k in snp_d:
        snp = snp_d[k][0]
        if snp.is_covered:
            snp.cleanse_counts(filtered_reads)
        
    pass

def show_data(transcript_d, snp_d):
    for k in sorted(transcript_d):
        transcript_shown = False
        # print transcript_d[k]
        show_snps = reduce(lambda x,y:x or y, [snp.is_covered for snp in transcript_d[k].snps])
        if not show_snps:
            continue        
        for snp in transcript_d[k].snps:
            if snp.is_covered:
                if not transcript_shown:
                    print transcript_d[k]
                    transcript_shown = True
                print '\t',
                print snp
#    for k in sorted(snp_d):
#        print k
#        for tx in snp_d[k][1:]:
#            print '\t' + tx
    pass    


def main(argv):
    
    transcript_d, snp_d = read_transcript_data(open(argv[0]))
    # show_data(transcript_d, snp_d)
    
    read_checklist = set([])
    for bam_fn in argv[1:]:
        process_pileups(pysam.Samfile(bam_fn, 'rb'), snp_d, read_checklist)    
        # show_data(transcript_d, snp_d)
        # break
    # print list(read_checklist), len(read_checklist)
    show_data(transcript_d, snp_d)
    
    pickle.dump(snp_d, open('SNPDATA.dat', 'wb'))
    pickle.dump(transcript_d, open('TRANSCRIPTDATA.dat', 'wb'))
    pickle.dump(read_checklist, open('READCHECKLIST.dat', 'wb'))
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

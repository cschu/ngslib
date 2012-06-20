#!/usr/bin/python

import os
import re
import sys
import math
import glob
import time

import pysam

from snptools import SNPLine, SNPDict


class IndexOutOfBoundsException(Exception):
    pass

#
class FASTASeq(object):
    def __init__(self, open_fn):
        self.id = open_fn.readline()[1:].strip()
        self.sequence = open_fn.read().replace('\n', '')
        pass
    def at(self, index, length=1, offset=0):  
        # TODO: needs range check for the whole sequence
        if index - offset >= len(self.sequence):
            errmsg = 'Pos %i (offset=%i) is outside sequence (length=%i)'
            raise IndexOutOfBoundsException(errmsg % (index, offset, 
                                                      len(self.sequence)))
        return self.sequence[index - offset: (index - offset) + length]
    pass

def median(lst):
    mid = len(lst)/2
    if len(lst) % 2 == 0:
        return sum(lst[mid-1:mid+1])/2.0
    else:
        return float(lst[mid])


def compare_sequences(positions, readseq, genome, cigar):
    ppos = 0
    spos = 0


    mismatches = []
    for op, length in cigar:
        if op == 0:
            cur_pos = [p for p in positions[ppos:ppos+length]]
            cur_seq = [c for c in readseq[spos:spos+length]]
            refseq = [genome.at(p) for p in cur_pos]
            
            cmpseq = zip(refseq, cur_seq)
            # print cmpseq
            miscount = sum([float(x[0]!=y[0]) 
                            for x, y in cmpseq])
            mismatches.append(miscount)
            spos += length
            ppos += length
        elif op == 1:
            spos += length
        elif op == 2:
            pass
        elif op == 3:
            pass
        elif op == 4:
            spos += length
        elif op == 5:
            spos += length
    print sum(mismatches)/len(readseq)
            


    pass


def parse_MD(mdstr):
    
    match_pos = 0
    mismatch_pos = 0
    while True:
        if not mdstr: break
        match = re.match('[0-9]+', mdstr)
        if match:
            match_pos += int(match.group())
            mdstr = mdstr[match.end():]
        else:
            match = re.match('([A-Z]|\^[A-Z]+)', mdstr)
            if match:
                mismatch_pos += len(match.group().strip('^'))
                mdstr = mdstr[match.end():]
            else:
                raise Exception('Weird mdstr: %s!' % mdstr)
        pass
    return mismatch_pos



def compute_seqerror(samfn, genome, snpdict, genomeid):
    samfile = pysam.Samfile(samfn, 'rb')
    seq_coverage = {}
        
    seqerrors = []
    for read in samfile.fetch(genomeid):
        
        """
        mismatches = []
        seqerr = compare_sequences(map(int, read.positions),
                                   read.seq,
                                   genome,
                                   read.cigar)
        print 'SEQERR', seqerr
        """                        
        """
        for pos, base in zip(map(int, read.positions), list(read.seq)):
            if base != genome.at(pos): 
                mismatches.append((pos, base, genome.at(pos)))
        seqerr = len(mismatches)/float(read.qlen)
        seqerrors.append(seqerr)
        # print mismatches, read.qlen, seqerr
        """
        aux_d = dict(read.tags)
        # seqerrors.append(aux_d['NM']/float(read.qlen))
        seqerr = parse_MD(aux_d['MD'])/float(read.qlen)

        seqerrors.append(seqerr)
        
        

        
        # print read.positions
        print samfn, genomeid
        print read.qname, read.qlen
        # print mismatches, len(mismatches), read.qlen, read.is_reverse
        print read.cigar
        # print read.qual
        # print read.qqual
        print read.tags
        print '***'

        if seqerr > 1.0:
            sys.stderr.write('WEIRD seqerr=%f. Exiting.\n')
            sys.exit(1)
        
        

    # sys.exit(1)
    # print 'SEQERRORS:', seqerrors, len(seqerrors), sum(seqerrors)
    return seqerrors



def main(argv):

    print 'START:', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    FASTAPATH = '/home/schudoma/tools/ngs/tair10'
    SNPFILE = '/home/schudoma/projects/ngs/filtered_variant.txt'
    snpdict = SNPDict(open(SNPFILE))    

    csv_out = open('seqerrors.dat', 'w')

    seqerr = []
    chrfiles = glob.glob(os.path.join(FASTAPATH, 'TAIR10_chr?.fas'))
    for chrfile in sorted(chrfiles):
        print 'Genome:', chrfile
        genome = FASTASeq(open(chrfile))
        genomeid = os.path.basename(chrfile).rstrip('.fas').lstrip('TAIR10_')
        genomeid = genomeid[0].upper() + genomeid[1:]
        for fn in glob.glob('*.concordant_uniq.bam'):            
            if fn.startswith('col-N-wot'): 
                continue
            print '\t%s' % fn
            print '\t%s' % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            
            seqerrors = compute_seqerror(fn, genome, snpdict, genomeid)
            seqerr.append((sum(seqerrors)/len(seqerrors),
                           median(sorted(seqerrors))))            

            for error in seqerrors:
                csv_out.write('%f\n' % error)

            print '\tSEQERROR=', seqerr[-1]
            # sys.exit(1)
            # seqerr.append((genomeid, fn, 
            #               compute_seqerror(fn, genome, snpdict, genomeid)))
            # print '\tSEQERROR=', seqerr[-1][2]

    # print 'AVERAGE SEQERR=', sum([x[2] for x in seqerr])/len(seqerr)
    csv_out.close()
    means = [x[0] for x in seqerr]
    medians = [x[1] for x in seqerr]
    print 'AVERAGE MEAN SEQERR=', sum(means)/len(means)            
    print 'AVERAGE MEDIAN SEQERR=', sum(medians)/len(medians)            

    print 'END:', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            
    return None
        


###
def compute_seqerror_per_pos(samfn, genome, snpdict, genomeid):
    samfile = pysam.Samfile(samfn, 'rb')
    seq_coverage = {}    
    i = 0
    for col in samfile.pileup(genomeid):

        counts = [read.alignment.seq[read.qpos] 
                  for read in col.pileups]
        counts = dict([(x, counts.count(x)) for x in set(counts)])        

        

        bases_at_pos = float(counts.get(genome.at(col.pos), 0))
        # print 'BASESATPOS_1:', bases_at_pos
        snp_query = snpdict.get((col.pos - 1, genomeid), None)
        if not snp_query is None:
            bases_at_pos += float(counts.get(snp_query.mutation, 0))      
        # print 'BASESATPOS_2:', bases_at_pos            
        
        if bases_at_pos/sum(counts.values()) > 1.0:
            print samfn, genomeid
            print genome.at(col.pos), col.pos
            print snp_query.mutation, snp_query.refbase
            print counts.items()
            sys.exit(1)
        seqerr = 1.0 - bases_at_pos/sum(counts.values())
        seq_coverage[col.pos] = seqerr

        # print 'coverage at base %s = %s' % (col.pos, col.n)
        # print chr1.at(col.pos)
        # print chr1.at(col.pos, offset=1)
        pass

        

    """
    fo = open('CHR1_err.csv', 'w')
    for pos, err in sorted(seq_coverage.items()):
        fo.write('%i,%f\n' % (pos, err))
    fo.close()
    """
    return sum(seq_coverage.values())/len(seq_coverage.values())




if __name__ == '__main__': main(sys.argv[1:])




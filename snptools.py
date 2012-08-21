#!/usr/bin/python

import os
import sys 
import math

#
genome_d = {'1': '1', 
            '2': '2',
            '3': '3',
            '4': '4',
            '5': '5',
            'C': 'C',
            'M': 'M',
            '6': 'C',
            '7': 'M'}


#
class SNPLine(object):
    """
    Modeled after "filtered_variant.txt @ 
    http://1001genomes.org/data/MPI/MPICao2010/releases/current/README
    Example line:
    Ped-0	1	110	G	T	27	3	1	1
    """
    def __init__(self, line, lineid=1):
        line = [col.strip() for col in line.split('\t')]
        self.seqid = 'Chr%c' % genome_d[line[1]]
        self.type_ = 'snp'
        self.pos = int(line[2])
        self.quality = int(line[5])
        try:
            self.supporting_reads = int(line[6])
        except:
            #print 'lineproblem:', line
            #sys.exit()
            self.supporting_reads = -1
            pass
        try:
            self.concordance = float(line[7])
        except:
            self.concordance = -1.0
        try:
            self.avg_align_per_overl_reads = float(line[8])
        except:
            self.avg_align_per_overl_reads = -1.0
        self.lineid = lineid
        self.refbase = line[3]
        self.mutation = line[4]
        self.source = line[0]
        pass
    def __repr__(self):
        outstring = [self.source, self.seqid[-1], self.pos,
                     self.refbase, self.mutation, self.quality,
                     self.supporting_reads, self.concordance,
                     self.avg_align_per_overl_reads]
                     
        return '\t'.join(map(str, outstring))
    def get_gffstring(self):        
        attributes = [
            ('quality', str(self.quality)), 
            ('supporting_reads', str(self.supporting_reads)),
            ('concordance', str(self.concordance)),
            ('avg_align_overl_reads', str(self.avg_align_per_overl_reads)),
            ('refbase', self.refbase),
            ('mutation', self.mutation)
            ]
        # 21-08-2012 was: ..., str(self.pos), str(self.pos), -> changed to UCSC format
        # according to http://code.google.com/p/bedtools/wiki/Usage#intersectBed
        fields = [self.seqid, '.', self.type_, 
                  str(self.pos - 1), str(self.pos), '.', '.', '.',
                  ';'.join(['%s=%s' % attr for attr in attributes])]
        return '\t'.join(fields)

#
class SNPDict(dict):
    def __init__(self, open_fn):        
        for line in open_fn:
            snpline = SNPLine(line, lineid=len(self) + 1)
            self[(snpline.seqid, snpline.pos)] = snpline
        pass
    pass




###
def main(argv):
    return None

if __name__ == '__main__': main(sys.argv[1:])

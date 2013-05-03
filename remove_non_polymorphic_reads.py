#!/usr/bin/env python
'''
Created on Apr 19, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

import pysam

from TranscriptIO import read_transcript_data, get_timestamp

def get_mapping_reads(bamfile, snp_d):
    N = len(snp_d)
    i = 1.0    
    valid_reads = set()
    
    for k in snp_d:
        if (i*100/N) % 5 == 0: 
            sys.stderr.write('%s: %i SNPs processed (%i%%)\n' % (get_timestamp(), i, i/N * 100))
        i += 1.0
        snp = snp_d[k][0]
        
        for read in bamfile.fetch(region=snp.get_region(k[0])):
            if read.flag & 0x4 == 0:
                valid_reads.add(read.qname)        
        
        return valid_reads
    
def write_mapping_reads(bamfile):
    pass

def main(argv):
    
    sys.stderr.write('%s: Reading transcript/SNP data\n' % get_timestamp())    
    transcript_d, snp_d = read_transcript_data(open(argv[0]))
    
    for bam_fn in argv[1:]:  
        sys.stderr.write('%s: Processing file %s...\n' % (get_timestamp(), bam_fn))      
        infile = pysam.Samfile(bam_fn, 'rb')
        mapping_reads = get_mapping_reads(infile, snp_d)    
        bam_fn_out = bam_fn.replace('.bam', '.prefiltered.bam')
        outfile = pysam.Samfile(bam_fn_out, 'wb', template=infile)
        
        for read in infile.fetch():
            if read.qname in mapping_reads:
                outfile.write(read)
        outfile.close()
        infile.close()
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python
'''
Created on Dec 6, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os

import analyse_snps_pe2 as analyse


def main(argv):
    version = '002'
    
    args = argv[1:]
    hit_mode = argv[0]
    if hit_mode not in ('inter', 'intra'):
        sys.stderr.write('Cannot find valid hit mode: %s. Exiting.\n' % hit_mode)
        sys.exit(1)
        
    all_snps, all_regions = {}, {}
    
    for bam_fn in args:
        snp_fn = bam_fn[:bam_fn.find('.bam')] + '_hits_snp_' + hit_mode + '.gff'
               
        sampleID = os.path.basename(bam_fn)
        sampleID = sampleID[:sampleID.find('.')] 
        print bam_fn, snp_fn, sampleID
        
        outfile = '%s.%s.%s.csv' % (version, sampleID, hit_mode)
        analysed_snps, reads_per_gene = analyse.analyse_snps(pysam.Samfile(bam_fn, 'rb'), 
                                                             open(snp_fn, 'rb'), 
                                                             hit_mode, 
                                                             out=open(outfile, 'ab'))
        for snp in analysed_snps:
            key = (snp[0], snp[1])
            all_snps[]
                
            
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

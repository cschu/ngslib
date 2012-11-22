#!/usr/bin/env python

'''
Created on Aug 23, 2012

@author: schudoma
'''

import sys
import subprocess

import gff_helpers

def main(argv):
    
    bamfile = argv[0]
    
    snp_gff = '/home/schudoma/projects/ngs/ped-0-snps-intergenic.gff'
    snps = gff_helpers.read_snp_from_gff(open(snp_gff))
    
    genome_path = '/home/schudoma/projects/ngs/tair10/TAIR10_chr%c.fas'
    
    base_cmd = ['samtools', 'mpileup', '-f']
    for snp in snps:
        # (gffline[0], int(gffline[3]) + 1, int(gffline[4]), comments['refbase'], comments['mutation'])
        
        genome_ref = genome_path % snp[0][-1]
        cmd = base_cmd + [genome_ref, '-r', '%s:%i-%i' % (snp[0], snp[1]-1, snp[1]+1), bamfile]        
        
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output = p.communicate()[0].strip()#.split('\n')
        if len(output) > 0:
            sys.stdout.write('%s:%i-%i:%c%c\n' % (snp[0], snp[1]-1, snp[1]+1, snp[3], snp[4]))
            sys.stdout.write('%s\n' % output)
            sys.stdout.flush()
        
        #if output[-1] != '<mpileup> Set max per-file depth to 8000':
        #    print output[2:]
        
    
    # samtools mpileup -gf ~/projects/ngs/tair10/TAIR10_chr1.fas -r Chr1:10264-10265 2012-08-22/ped-N-shoot.all.sorted.bam
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
    
    
#!/usr/bin/env python
'''
Created on Jan 17, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import csv

from analyse_snps_pe3 import compute_consensus


def main(argv):
    
    out = sys.stdout
    col_primers, ped_primers = {}, {}
    
    for fn in argv:
        #sys.stdout.write('Processing' + fn + '\n')
        #sys.stdout.flush()
        reader = csv.reader(open(fn, 'rb'), delimiter=',', quotechar='"')
        for row in reader:
            if row[0] not in ('AGI', 'Contig'):
                key = (row[0], row[1])
                col, ped = row[12:14]
                if col != '':
                    col_primers[key] = col_primers.get(key, []) + [col]
                if ped != '':
                    ped_primers[key] = ped_primers.get(key, []) + [ped]
        pass
    
    #col_primers = comput(col_primers)
    #ped_primers = merge_primers(ped_primers)
    
    out.write('Contig/AGI,Position,Col-SNP,Ped-SNP\n')
    for key in sorted(set(col_primers.keys() + ped_primers.keys())):
        colseq = compute_consensus(col_primers.get(key, []))
        pedseq = compute_consensus(ped_primers.get(key, []))
        
        out.write('%s,%s,%s,%s\n' % (key + (colseq, pedseq))) 
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

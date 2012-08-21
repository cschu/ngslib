#!/usr/bin/env python


'''
Created on Aug 21, 2012

@author: Chris
'''

import sys

from check_snps import read_snp_data

def main(argv):
    snp_data = read_snp_data(open(argv[0]))
    
    out = sys.stdout
    
        
    for snpline in sorted(snp_data.items()):
        snpline = snpline[1]
        if snpline.refbase in 'ACGUT' and snpline.mutation in 'ACGUT':                
            out.write('%s\n' % snpline.get_gffstring())
        else:
            # out.write('%s\n' % snpline.get_gffstring())
            sys.stderr.write('???: %s\n' % str(snpline))

    
    
    pass


if __name__ == '__main__': main(sys.argv[1:])


# a number of snplines (43?) have ambiguity codes 
#$ /cygdrive/c/Users/Chris/workspace/ngslib/convert_1001.py ped-0-snps_no-indels
#.txt > ped-0-snps_no-indels.gff
#???: Ped-0      1       7056151 A       N       1       1       1.0     1.0
#???: Ped-0      1       12192624        Y       C       40      31      1.01.0
#???: Ped-0      1       12294696        T       N       2       2       2.01.0
#???: Ped-0      1       14196292        C       N       1       1       1.01.0
#???: Ped-0      1       14196293        T       N       1       1       1.01.0
#???: Ped-0      1       15083672        T       N       0       1       inf2.0
#???: Ped-0      1       15084619        A       N       0       1       inf3.125
#
#???: Ped-0      1       16431565        A       N       1       1       1.01.0
#???: Ped-0      1       17111141        W       T       17      3       1.02.75
#???: Ped-0      2       7521    G       N       3       4       1.33333 1.0
#???: Ped-0      2       85781   Y       T       32      27      1.0     1.0
#???: Ped-0      2       85834   K       T       32      22      0.9565221.0
#???: Ped-0      2       85836   S       G       40      22      1.0     1.0
#???: Ped-0      2       85844   K       T       40      25      1.0     1.0
#???: Ped-0      2       86536   K       T       26      4       1.0     1.0
#???: Ped-0      2       3627608 T       N       0       1       inf     3.6
#???: Ped-0      2       7492720 N       G       40      14      1.0     1.0
#???: Ped-0      2       11905490        K       T       40      10      1.01.0
#???: Ped-0      2       11905532        M       C       40      11      1.01.0
#???: Ped-0      2       11905533        S       G       40      10      1.01.0
#???: Ped-0      2       11905548        K       G       32      7       1.01.0
#???: Ped-0      2       11905549        R       A       32      7       1.01.0
#???: Ped-0      2       11905580        R       A       40      9       1.01.0
#???: Ped-0      2       11905582        Y       C       40      9       1.01.0
#???: Ped-0      2       11905770        R       A       36      3       1.01.0
#???: Ped-0      2       11905771        M       C       36      4       1.01.0
#???: Ped-0      2       11905772        Y       T       36      4       1.01.0
#???: Ped-0      2       11905773        W       A       36      4       1.01.0
#???: Ped-0      2       11905850        R       A       34      6       1.01.0
#???: Ped-0      2       11905852        W       A       36      5       1.01.0
#???: Ped-0      2       11905853        W       T       38      5       1.01.0
#???: Ped-0      2       11905854        W       A       30      5       0.833333
#        1.0
#???: Ped-0      2       14155545        R       G       14      2       1.01.0
#???: Ped-0      2       14155547        W       A       13      2       1.01.0
#???: Ped-0      2       16303213        K       T       30      3       1.01.0
#???: Ped-0      2       19677649        M       C       40      39      1.01.0
#???: Ped-0      3       15128557        M       A       24      16      1.02.476
#19
#???: Ped-0      3       15128576        M       C       40      11      1.01.0
#???: Ped-0      3       15128612        Y       C       40      10      1.01.153
#85
#???: Ped-0      4       2960867 C       N       1       1       1.0     1.33333
#???: Ped-0      4       3950841 T       N       1       1       1.0     3.33333
#???: Ped-0      4       3951153 A       N       0       1       inf     4.42857
#???: Ped-0      4       3954000 A       N       0       1       1.0     3.85714
#???: Ped-0      4       12493262        M       A       34      4       1.01.0
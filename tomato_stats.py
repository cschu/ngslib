#!/usr/bin/env python

'''
Created on Aug 29, 2012

@author: schudoma
'''

import sys

def percentage_str(n, total):
    return '(%.3f%%)' % (float(n)/total * 100)

def process_set(setname, values):
    set_d = dict(values)
    
    total_reads = sum(set_d.values())
    
    nomapping = sum([v for k, v in set_d.items() if 'nomapping' in k])
    up_unique = sum([v for k,v in set_d.items() if 'unpaired_uniq' in k])
    p_unique = sum([v for k,v in set_d.items() if 'concordant_uniq' in k or '.paired_uniq_' in k])
    # print unique
    #up_unique = [item for item in unique if 'unpaired' in item[0]]
    #print up_unique 
    transloc = sum([v for k,v in set_d.items() if 'transloc' in k])
    up_mult = sum([v for k,v in set_d.items() if 'unpaired_mult' in k])
    p_mult = sum([v for k,v in set_d.items() if 'concordant_mult' in k or '.paired_mult' in k])
    
    h_unique = sum([v for k,v in set_d.items() if 'halfmapping_uniq' in k])
    h_mult = sum([v for k,v in set_d.items() if 'halfmapping_mult' in k])
    
    print setname
    print 'Total reads (after quality trimming):', total_reads
    print 'Reads without match to Tomato genome (ITAG 2.3):', nomapping, percentage_str(nomapping, total_reads)
    print 'Unpaired reads with unique match:', up_unique, percentage_str(up_unique, total_reads)
    print 'Paired reads with unique match:', p_unique, percentage_str(p_unique, total_reads)
    print 'Unpaired reads with multiple matches:', up_mult, percentage_str(up_mult, total_reads)
    print 'Paired reads with multiple matches:', p_mult, percentage_str(p_mult, total_reads)
    
    print 'Paired reads with only one partner matching (unique):', h_unique, percentage_str(h_unique, total_reads)
    print 'Paired reads with only one partner matching (multiple):', h_mult, percentage_str(h_mult, total_reads)
    print 'Reads with translocation:', transloc, percentage_str(transloc, total_reads)
    
    
    print nomapping+up_unique+up_mult+p_unique+p_mult+transloc+h_unique+h_mult == total_reads
    pass


def main(argv):
    
    data = [(line.split(',')[0], int(line.split(',')[1].strip()))             
            for line in open(argv[0]).readlines()]
    data=sorted(data)
    data_d = {}
    for k, value in data:
        key = '_'.join(k.split('_')[:2])
        data_d[key] = data_d.get(key, []) + [(k, value)]
    for k, v in data_d.items():
        process_set(k, v)
    
    
    
    pass    

if __name__ == '__main__': main(sys.argv[1:])

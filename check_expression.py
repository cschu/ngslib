#!/usr/bin/python

'''
Created on Aug 6, 2012

@author: schudoma
'''
import re
import sys
import string

from get_expression_data import format_agi

def main(argv):
    
    query_id_file = 'tair_ids.txt'
    expr_data_file = 'expression_data.csv'
    
    query_ids = map(format_agi, map(string.strip, open(query_id_file).readlines()))
    # query_ids = [qid.rstrip('.1') for qid in query_ids]
    query_ids = list(set([re.sub('\.[0-9]', '', qid) for qid in query_ids]))
    expr_data = open(expr_data_file).readlines()
    # expr_ids  = re.findall('At[1-5]g[0-9]{5}', expr_data)
    expr_ids = ','.join([line.split()[0] 
                         for line in expr_data 
                         if line.startswith('At')]).replace('|', ',').split(',')
    fo = open('EXPRDATA.txt', 'w')
    for eid in expr_ids:
        fo.write('%s\n' % eid)
    fo.close()
    
    # Test
    # print 'At5g60548' in expr_ids, 'At5g60550' in expr_ids
    # return None
    expr_ids_set = set(expr_ids)
    not_in_set = 0
    in_set = 0
    for qid in query_ids:
        if qid not in expr_ids_set:
            not_in_set += 1
            print qid
        else:
            in_set += 1
    
    sys.stderr.write('%i %i %i\n' % (len(expr_ids), len(expr_ids_set), len(query_ids)))
    sys.stderr.write('%i %i\n' % (in_set, not_in_set))
    sys.stderr.write('%i %i\n' % (len(query_ids), len(set(query_ids))))
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])




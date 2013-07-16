#!/usr/bin/env python
'''
Created on Jun 28, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import re
import sys

from collections import defaultdict, Counter


def get_opt_fields(fields):
    opt_fields = {}
    for field in fields:
        key, value = field[:field.find(':')], field[field.rfind(':') + 1:]
        opt_fields[key] = value
    return opt_fields

def main(argv):
    
    #matches = defaultdict(dict(fwd=list(), rev=list()))
    matches = defaultdict(defaultdict)
    for fn in argv:
        sys.stderr.write('Processing file %s...\n' % fn)
        c = 0
        for line in open(fn).readlines():
            if line.startswith('@'):
                continue
            # print c, line
            fields = line.strip().split()            
            flag = int(fields[1])
            
            if flag & 4 == 4:
                continue
            if re.search('NM:i:[012]\s', line):
                key = (fields[0], fields[2])
                
                if (flag & 64 == 64) or (flag & 1 == 0):
                    matches[key][1] = fields
                else:
                    matches[key][2] = fields
            #c += 1
            #if c == 10: break
    
    #for item in matches.items()[:10]:
    #    print item        
    #print matches.keys()[0], matches[matches.keys()[0]]['fwd']
    n_weird_orientation = 0    
    transcripts = defaultdict(defaultdict)
    sys.stderr.write('Checking transcripts...\n')    
    for key in matches.keys():
        if 1 in matches[key] and 2 in matches[key]:
            # both mates aligned, need to check
            flag1, flag2 = int(matches[key][1][1]), int(matches[key][2][1])
            if not((flag1 & 16 == 0 and flag2 & 16 == 16) or (flag1 & 16 == 16 and flag2 & 16 == 0)):
                # both mates need to have complementary orientation, else discard pair
                n_weird_orientation += 1
                del matches[key]
                continue
            optfields1, optfields2 = get_opt_fields(matches[key][1]), get_opt_fields(matches[key][2])
            # n_mismatches = 
            transcripts[key[1]]['pairs'] = transcripts[key[1]].get('pairs', []) + [(int(optfields1['NM']),
                                                                                    int(optfields2['NM']),)]
        else:
            mate = 1 if 1 in matches[key] else 2
            optfields = get_opt_fields(matches[key][mate])
            transcripts[key[1]]['singles'] = transcripts[key[1]].get('singles', []) + [int(optfields['NM'])]
    sys.stderr.write('%i pairs with questionable orientation found...\n' % n_weird_orientation)
    
    
    headers = ['AGI', '#singles', '#pairs'] + ['#singles_%i' for i in xrange(3)] + ['#pairs_%i' for i in xrange(5)]
    sys.stdout.write(','.join(headers) + '\n')
    for key in sorted(transcripts):
        pairs = Counter(map(lambda x:sum(x), transcripts[key]['pairs'])) if 'pairs' in transcripts[key] else Counter() 
        singles = Counter(transcripts[key]['singles']) if 'singles' in transcripts[key] else Counter()        
        n_pairs, n_singles = sum(pairs.values()), sum(singles.values())
        out_str = ','.join(map(str, [key, n_singles, n_pairs] + [singles[i] for i in xrange(3)] + [pairs[i] for i in xrange(5)]))
        sys.stdout.write(out_str + '\n')
    
    
    pass
    

if __name__ == '__main__': main(sys.argv[1:])

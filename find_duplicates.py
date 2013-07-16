#!/usr/bin/env python
'''
Created on Jun 20, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import re
import sys

from collections import defaultdict


def find_final_target(map_d, key):
    while True:
        # maps_to[sid1] = maps_to[sid2]
        if key not in map_d: break
        key = map_d[key]
    return key

def main(argv):
    
    maps_to = {}
    seq_d = {}
    
    fn = argv[0]
    for line in open(fn):
        if line.startswith('@'):
            continue
        fields = line.strip().split()
        is_mapped = int(fields[1]) & 4 == 0 
                 
        if not is_mapped:
            continue
        
        sid1 = fields[0]
        sid2 = fields[2].strip('@')
        seq_d[sid1] = fields[9]
        
        if sid1 in maps_to:
            # must not happen, due to bowtie-mode!
            sys.stderr.write("WEIRD: %s already known\n" % sid1)
            sys.exit(1)        
        
        allowed_match = re.search('NM:i:0', line)
        if allowed_match:
            if sid2 in maps_to:
                # if sid1 -> sid2 and sid2 -> sidX, then sid1 -> sidX
                maps_to[sid1] = find_final_target(maps_to, sid2)
            else:
                # sid1 -> sid2
                maps_to[sid1] = find_final_target(maps_to, sid2)
                for sid in maps_to:
                    if maps_to[sid] == sid1:
                        maps_to[sid] = maps_to[sid1]
            #if sid2 in maps_to:
            #    # if sid1 -> sid2 and sid2 -> sidX, then sid1 -> sidX
            #    maps_to[sid1] = maps_to[sid2]
            #else:
            #    # sid1 -> sid2
            #    maps_to[sid1] = sid2
        else:
            maps_to[sid1] = sid1
            pass 
    
    mapped_by = defaultdict(set)
    for k in maps_to:
        mapped_by[maps_to[k]].add(k)
    
    dummyset = set()
    for k in mapped_by:
        mapped_by[k] = list(set(mapped_by[k]) - set([k]))
        print '>', k, len(mapped_by[k]) + 1
        seq = sorted([(len(seq_d[kk]), seq_d[kk]) for kk in [k] + mapped_by[k] if kk in seq_d])[-1]
        print 'SEQ:', seq[1]
        for kk in mapped_by[k]:
            print '\t' + kk
            dummyset.add(kk)
        dummyset.add(k)
    sys.stderr.write('DUMMYSET: %i ELEMENTS\n' % len(dummyset))
    
    
    
    pass



def main2(argv):
    
    fn = argv[0]
    
    adj_list = defaultdict(set)
    
    for line in open(fn):
        if line.startswith('@'):
            continue
        fields = line.strip().split()
        is_mapped = int(fields[1]) & 4 == 0 and re.search('NM:i:0', line)

        sid1 = fields[0]        
        if is_mapped:            
            sid2 = fields[2].strip('@')
        
            if sid1 == sid2:
                adj_list[sid1] = adj_list[sid1].union(set())
            else:
                adj_list[sid1].add(sid2)
                adj_list[sid2].add(sid1)
        else:
            adj_list[sid1] = adj_list[sid1].union(set())
    
    
    out = open('SEQUENCE_CLUSTERS.txt', 'wb')
    while len(adj_list) > 0:
        max_deg = max([len(adj_list[sid]) for sid in adj_list])
        sid_max = sorted([sid for sid in adj_list if len(adj_list[sid]) == max_deg])[0]
        sid_max_adj = set().union(adj_list[sid_max])
        print sid_max, len(sid_max_adj)
        print sid_max_adj
        
        del adj_list[sid_max]
        for sid in sid_max_adj:
            print '**', sid, len(sid_max_adj)
            
            print '##', adj_list[sid]
            adj_list[sid] = adj_list[sid].difference(sid_max_adj)            
            print '###', adj_list[sid]
            
            if len(adj_list[sid]) == 0:
                del adj_list[sid]
                print '***', sid, len(sid_max_adj)
            else:
                sid_max_adj = sid_max_adj.difference(set([sid]))
                print '****', sid, len(sid_max_adj)
            sys.exit(1)
        
        out.write('>%s %i\n' % (sid_max, len(sid_max_adj) + 1))
        for sid in sid_max_adj:
            out.write(' %s\n' % sid)
        out.flush()
        pass
    out.close()
        
                
        
        
        
        
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

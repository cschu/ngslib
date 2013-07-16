#!/usr/bin/env python
'''
Created on Jul 11, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import time, datetime

def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')


def main(argv):
    
    sys.stderr.write('%s> Generating set1 from file %s...\n' % (get_timestamp(), argv[0]))
    set1 = set([line.strip() for line in open(argv[0]) if line[0] != '*'])
    sys.stderr.write('%s> Generating set2 from file %s...\n' % (get_timestamp(), argv[1]))
    set2 = set([line.strip() for line in open(argv[1]) if line[0] != '*'])
    out = open('%s_MOBILE_CANDIDATES.txt' % argv[2], 'wb')
    # out = sys.stdout    
    
    sys.stderr.write('%s> Intersecting sets...\n' % get_timestamp())
    # intersection = sorted(set1.intersection(set2))
    intersection = list(set1.intersection(set2))
    sys.stderr.write('%s> Sorting intersection...\n' % get_timestamp())
    intersection.sort()
    
    last_id = None
    sets = []    
    sys.stderr.write('%s> Writing fq-ids...\n' % get_timestamp())
    for id_ in intersection:
        id_, set_id = id_.split('_')
        if last_id is None or id_ == last_id:
            pass
        else:
            out.write('%s\t%s\n' % (last_id, ','.join(sorted(sets))))
            sets = []
        sets.append(set_id)
        last_id = id_
    out.write('%s\t%s\n' % (last_id, ','.join(sorted(sets))))
    sys.stderr.write('%s> Done.\n' % get_timestamp())
    out.close()
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

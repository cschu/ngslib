#!/usr/bin/env python
'''
Created on Jul 25, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re
import os
import time, datetime

def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
    
is_match = re.compile('NM:i:[01]\s')

def main(argv):
    
    mapping, not_mapping = set(), set()
    hostname = argv[0]
    output = argv[1]
    
    for fn in argv[2:]:
        sys.stderr.write('%s> Processing file %s...\n' % (get_timestamp(), fn))
        setid = os.path.basename(fn).split('_')[1]
        sampleid = os.path.basename(fn).split('_')[0]        
        sys.stderr.write('%s> Reading...\n' % get_timestamp()) 
        for line in open(fn):
            if line.startswith('@'):
                continue
            fields = line.strip().split()
            is_mapped = (int(fields[1]) & 4 == 0) and is_match.search(line)
            read_id = fields[0] + '_' + setid
            
            if is_mapped:
                not_mapping.discard(read_id)
                mapping.add(read_id)
            elif read_id not in mapping:
                not_mapping.add(read_id)
            else:
                pass
            
            
    write_mapping = output == 'MAPPING' or output == 'BOTH'
    write_notmapping =  output == 'NOTMAPPING' or output == 'BOTH'
            
    if write_mapping:            
        sys.stderr.write('%s> Writing mapped reads...\n' % (get_timestamp()))
        fn_out = '%s_VS_%s_MAPPING.txt' % (sampleid, hostname)
        out = open(fn_out, 'wb')
        for mate_id in sorted(mapping):
            out.write('%s\n' % mate_id)
        out.close()
    if write_notmapping:
        sys.stderr.write('%s> Writing unmapped reads...\n' % (get_timestamp()))
        fn_out = '%s_VS_%s_NOTMAPPING.txt' % (sampleid, hostname)
        out = open(fn_out, 'wb')
        for mate_id in sorted(not_mapping):
            out.write('%s\n' % mate_id)
        out.close()
    sys.stderr.write('%s> Done.\n' % (get_timestamp()))
        
                
               
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

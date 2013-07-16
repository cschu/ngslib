#!/usr/bin/env python
'''
Created on Jun 18, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import glob
import re
import os

import time, datetime

#from collections import defaultdict

"""
ACL1_ON_CCL1/R1_paired_on_pR1_bt2index.sam  ACL1_ON_CCL1/unpaired_all_on_pR1_bt2index.sam  CCL1_ON_ACL1/R2_paired_on_pR1_bt2index.sam
ACL1_ON_CCL1/R1_paired_on_pR2_bt2index.sam  ACL1_ON_CCL1/unpaired_all_on_pR2_bt2index.sam  CCL1_ON_ACL1/R2_paired_on_pR2_bt2index.sam
ACL1_ON_CCL1/R1_paired_on_u_bt2index.sam    ACL1_ON_CCL1/unpaired_all_on_u_bt2index.sam    CCL1_ON_ACL1/R2_paired_on_u_bt2index.sam
ACL1_ON_CCL1/R2_paired_on_pR1_bt2index.sam  CCL1_ON_ACL1/R1_paired_on_pR1_bt2index.sam     CCL1_ON_ACL1/unpaired_all_on_pR1_bt2index.sam
ACL1_ON_CCL1/R2_paired_on_pR2_bt2index.sam  CCL1_ON_ACL1/R1_paired_on_pR2_bt2index.sam     CCL1_ON_ACL1/unpaired_all_on_pR2_bt2index.sam
ACL1_ON_CCL1/R2_paired_on_u_bt2index.sam    CCL1_ON_ACL1/R1_paired_on_u_bt2index.sam       CCL1_ON_ACL1/unpaired_all_on_u_bt2index.sam
"""
def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')


def main(argv):
    
    mapping = set([])
    not_mapping = set([])    
    
    # we're interested in mobile RNAs, hence we only check for
    # reads from the parasite-set (i.e., in e.g. CAL vs CM, 
    # CAL is assumed to be the parasite!) -- ATC must always be host!!!
    parasite = argv[0]
    host = argv[1]
    parasite_bn = os.path.basename(parasite).split('_')[0]
    host_bn = os.path.basename(host).split('_')[0]    
    hostfiles = glob.glob(os.path.join(host, '*.sam'))
    parasitefiles = glob.glob(os.path.join(parasite, '*.sam'))
    
    mismatch_pattern = re.compile('NM:i:[01]\s')
    
    # Since we're mapping both forward (i.e. parasite vs host) and backward (i.e. host vs. parasite),
    # it is possible that a read is marked both as non-mapping and as mapping.
    # We do two-step-mapping because long query reads will not align against reads that would match but 
    # are shorter than the query (e.g. due to clipping).    
    for fn in (hostfiles + parasitefiles):
        sys.stderr.write('%s> Processing file %s...\n' % (get_timestamp(), fn))
        is_parasitefile = fn in parasitefiles
        if 'R1_paired' in fn:
            id_suffix1 = '_R1'
        elif 'R2_paired' in fn:
            id_suffix1 = '_R2'
        else:
            id_suffix1 = '_U'
        if 'pR1' in fn:
            id_suffix2 = '_R1'
        elif 'pR2' in fn:  
            id_suffix2 = '_R2'
        else:
            id_suffix2 = '_U'
        
        sys.stderr.write('%s> Reading...\n' % get_timestamp())        
        for line in open(fn):
            if line.startswith('@'):
                continue
            
            fields = line.strip().split()
            is_mapped = (int(fields[1]) & 4) == 0 and mismatch_pattern.search(line) is not None
                    
            if is_parasitefile:
                mate_id = fields[0].strip() + id_suffix1
            else:
                mate_id = fields[2].strip().strip('@') + id_suffix2
        
            if is_mapped:
                not_mapping.discard(mate_id)
                mapping.add(mate_id)
            elif mate_id not in mapping:
                not_mapping.add(mate_id)
            else:
                pass
            
    # this is no longer needed
    # sys.stderr.write('%s> Set sizes: mapping:%i not_mapping:%i\n' % (get_timestamp(), len(mapping), len(not_mapping)))      
    # sys.stderr.write('%s> Adjusting sets...\n' % (get_timestamp()))
    # not_mapping.difference_update(mapping)
    # sys.stderr.write('%s> Set sizes: mapping:%i not_mapping:%i\n' % (get_timestamp(), len(mapping), len(not_mapping)))
    
    sys.stderr.write('%s> Writing mapped reads...\n' % (get_timestamp()))
    fn_out = '%s_VS_%s_MAPPING.txt' % (parasite_bn, host_bn)
    out = open(fn_out, 'wb')
    for mate_id in sorted(mapping):
        out.write('%s\n' % mate_id)
    out.close()
    
    sys.stderr.write('%s> Writing unmapped reads...\n' % (get_timestamp()))
    fn_out = '%s_VS_%s_NOTMAPPING.txt' % (parasite_bn, host_bn)
    out = open(fn_out, 'wb')
    for mate_id in sorted(not_mapping):
        out.write('%s\n' % mate_id)
    out.close()
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python
'''
Created on Jun 20, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os

from collections import defaultdict

"""
Sample_CCL1/CCL1_SPECIFIC.txt  Sample_CCL2/CCL2_SPECIFIC.txt <- refsets
CCL1_VS_ACL1_MAPPING.txt  CCL1_VS_ACL2_MAPPING.txt <- checksets

HWI-ST377:249:C1JUVACXX:8:1101:10000:100455_R1
HWI-ST377:249:C1JUVACXX:8:1101:10000:100455_R2
HWI-ST377:249:C1JUVACXX:8:1101:10000:14535_R1
HWI-ST377:249:C1JUVACXX:8:1101:10000:14535_R2
HWI-ST377:249:C1JUVACXX:8:1101:10000:16421_R1
HWI-ST377:249:C1JUVACXX:8:1101:10000:16421_R2
HWI-ST377:249:C1JUVACXX:8:1101:10000:18304_R1
HWI-ST377:249:C1JUVACXX:8:1101:10000:18304_R2
HWI-ST377:249:C1JUVACXX:8:1101:10000:19530_R1
HWI-ST377:249:C1JUVACXX:8:1101:10000:19530_R2
 
HWI-ST377:249:C1JUVACXX:8:1101:10000:100455
HWI-ST377:249:C1JUVACXX:8:1101:10000:12265
HWI-ST377:249:C1JUVACXX:8:1101:10000:14535
HWI-ST377:249:C1JUVACXX:8:1101:10000:16421
HWI-ST377:249:C1JUVACXX:8:1101:10000:18304
HWI-ST377:249:C1JUVACXX:8:1101:10000:19530
HWI-ST377:249:C1JUVACXX:8:1101:10000:31570
HWI-ST377:249:C1JUVACXX:8:1101:10000:33222
HWI-ST377:249:C1JUVACXX:8:1101:10000:36902
HWI-ST377:249:C1JUVACXX:8:1101:10000:44719
"""
#Sample_CAL1/CAL1_ATH_CONSERVED.txt CCL1_VS_CAL1_NOTMAPPING.txt

def main(argv):
    
    fn_refset = argv[0]
    fn_checkset = argv[1]
      
    #fn_out = '%s_IN_%s_MOBILE_CANDIDATES.txt' % (os.path.basename(fn_refset)[:4],)
    fn_out = '%s_MOBI_CAND.txt' % os.path.basename(fn_checkset).replace('VS', 'IN')
    
    refset = set()
    found_ids = defaultdict(list)
    
    for line in open(fn_refset):
        refset.add(line.strip())
    
    out = open(fn_out, 'wb')
    for line in open(fn_checkset):
        fqid, mid = line.strip().split('_')
        if fqid in refset:
            # out.write('%s\n' % line.strip())
            if mid not in ('R1', 'R2'):
                found_ids[fqid] = []
            else:
                found_ids[fqid].append(mid)
    
    for fqid in sorted(found_ids):
        #out.write('%s\n' % line.strip())
        out.write(' '.join([fqid] + sorted(found_ids[fqid])) + '\n')
        out.flush()
    out.close()
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

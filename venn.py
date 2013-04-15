#!/usr/bin/env python
'''
Created on Apr 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import os
import sys
import cPickle

from Transcript import Transcript
from SNP_Position import SNP_Position
from TranscriptIO import load_transcript_data

"""
Dickes B, home an der Spree,
im Sommer tust du gut und im Winter tut's weh.
Mama Berlin - Backsteine und Benzin
- wir lieben deinen Duft, wenn wir um die Haeuser ziehn.

Dickes G, paarst dich mit dem C
manchmal auch mit U aber niemals mit G
Nur wenn auch Hoogsteen oder Sugar mit ziehn
 - dann gibt's auch mal ein (Paar aus) Homo-Guanin.
"""



def vennAB(A, B):    
    return [(('%s') % A[0], A[1] - B[1]), 
            (('%s') % A[0], B[1] - A[1]), 
            ('%s_+_%s' % (A[0], B[0]), A[1].intersection(B[1]))]

def vennABC(A, B, C):    
    return [(('%s') % A[0], A[1] - C[1].union(B[1])), 
            (('%s') % B[0], B[1] - A[1].union(C[1])), 
            (('%s') % C[0], C[1] - B[1].union(A[1])),
            ('%s_+_%s' % (A[0], B[0]), A[1].intersection(B[1]) - C[1]),
            ('%s_+_%s' % (A[0], C[0]), A[1].intersection(C[1]) - B[1]),
            ('%s_+_%s' % (B[0], C[0]), B[1].intersection(C[1]) - A[1]),
            ('%s_+_%s+%s' % (A[0], B[0], C[0]), A[1].intersection(B[1]).intersection(C[1]))]

def load_set(data, use_all=True):
    set_ = set()
    for k in data:
        if data[k].has_covered_snps():            
            if use_all or (data[k].is_mobile and data[k].binom_score > 10.0 ** -5):
                set_.add(k)
    return set_

def process_sets(sets):
    if len(sets) == 2:
        return vennAB(sets[0], sets[1])
    elif len(sets) == 3:
        return vennABC(sets[0], sets[1], sets[2])
    else:
        return {}
    
        
"""
def process_sets(sets):
    d = {i: set_i[1] for i, set_i in enumerate(sets)}
    for i, set_i in enumerate(sets):
        for j in d:
            if i != j and i not in j:
                
    
    
    
    for i, set_i in enumerate(sets):
        for j, set_j in enumerate(sets):
            if j > i:
                
    
    pass
"""

def main(argv):
    sets = []
    use_mobile = argv[0] == '--USE_MOBILE'
    for fn in argv[1:]:
        if 'ped' in os.path.basename(fn).lower():
            sample = 'Ped'
        else:
            sample = 'Col'
        
        set_ = load_set(load_transcript_data(open(fn, 'rb'), sample=sample), not use_mobile)
        sets.append((os.path.basename(fn).replace('_TRANSCRIPTDATA.pickled', ''), set_))
    
    print 'SETS =', len(sets)
    print ' '.join(['%s: %i' % (s[0], len(s[1])) for s in sets])
    
    for id_, set_ in process_sets(sets):
        print id_, len(set_)
    
    for id_, set_ in process_sets(sets):
        print id_, len(set_)
        for item in set_:
            print item
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

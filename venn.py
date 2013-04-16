#!/usr/bin/env python
'''
Created on Apr 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import os
import sys
import cPickle

from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
from matplotlib import pyplot as plt

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
            (('%s') % B[0], B[1] - A[1]), 
            ('%s_+_%s' % (A[0], B[0]), A[1].intersection(B[1]))]

def vennABC(A, B, C):    
    return [(('%s') % A[0], A[1] - C[1].union(B[1])), 
            (('%s') % B[0], B[1] - A[1].union(C[1])), 
            (('%s') % C[0], C[1] - B[1].union(A[1])),
            ('%s_+_%s' % (A[0], B[0]), A[1].intersection(B[1]) - C[1]),
            ('%s_+_%s' % (A[0], C[0]), A[1].intersection(C[1]) - B[1]),
            ('%s_+_%s' % (B[0], C[0]), B[1].intersection(C[1]) - A[1]),
            ('%s_+_%s_+_%s' % (A[0], B[0], C[0]), A[1].intersection(B[1]).intersection(C[1]))]
    
def vennN(sets):
    # V = [(set([id_]), set_) for id_, set_ in sets]
    V = {i: (set([id_]), set_) for i, (id_, set_) in enumerate(sets)}
    setnames = {(id_,) : i for i, (id_, set_) in enumerate(sets)}
        
    print V
    p = 0
    while True:
        if p >= len(V): break        
        id1, set1 = V[p]
        # print 'p=', p, id1, set1        
        U = []
        for q in sorted(V):
            (id2, set2) = V[q]            
            if q > p and len(set(id1).intersection(set(id2))) == 0:
                # print 'q=', q, id2, set2            
                new_id = tuple(sorted(id1.union(id2)))
                new_set = set1.intersection(set2)                
                if new_id in setnames:
                    # can this occur without the existing set being empty?
                    # print 'Old set:', new_id
                    new_set = new_set.union(V[setnames[new_id]][1])
                    V[setnames[new_id]][1].update(new_set)                    
                else:                    
                    # print 'New set:', new_id
                    setnames[new_id] = len(V) + len(U)
                    U.append((len(V) + len(U), (new_id, new_set)))
                                        
                set1.difference_update(new_set)
                set2.difference_update(new_set)               
                # print 'U=', U
        V.update(dict(U))
        # print 'V=', V
        # print
        p += 1
    return V
        
                    
        
    

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
return [(('%s') % A[0], A[1] - C[1].union(B[1])), 
            (('%s') % B[0], B[1] - A[1].union(C[1])), 
            (('%s') % C[0], C[1] - B[1].union(A[1])),
            ('%s_+_%s' % (A[0], B[0]), A[1].intersection(B[1]) - C[1]),
            ('%s_+_%s' % (A[0], C[0]), A[1].intersection(C[1]) - B[1]),
            ('%s_+_%s' % (B[0], C[0]), B[1].intersection(C[1]) - A[1]),
            ('%s_+_%s_+_%s' % (A[0], B[0], C[0]), A[1].intersection(B[1]).intersection(C[1]))]
"""

def generate_diagram(sets, fn):
    fig = plt.figure(figsize=(4, 4))
    if len(sets) == 3:
        diagram_f = venn2
        circles_f = venn2_circles      
        subsets = [len(sets[0][1]), len(sets[1][1]), len(sets[2][1])]
        set_labels = map(lambda x: x.replace('Sample_', ''), [sets[0][0], sets[1][0]])
    elif len(sets) == 7:
        diagram_f = venn3
        circles_f = venn3_circles
        subsets = [len(sets[0][1]), len(sets[1][1]), len(sets[3][1]), len(sets[2][1]), 
                   len(sets[4][1]), len(sets[5][1]), len(sets[6][1])]
        set_labels = map(lambda x: x.replace('Sample_', ''), [sets[0][0], sets[1][0], sets[2][0]])
    else:
        pass 
    
    v = diagram_f(subsets=subsets, set_labels=set_labels)
    c = circles_f(subsets=subsets, linestyle='dashed')
    # plt.title("Sample Venn Diagram")    
    fig.savefig(fn)
    pass


def main(argv):
    sets = []
    use_all = True
    if argv[0] == '--USE_MOBILE':
        use_all = False
        argv = argv[1:]
    
    setnames = []
    for fn in argv:
        if 'ped' in os.path.basename(fn).lower():
            sample = 'Ped'
        else:
            sample = 'Col'
        
        set_ = load_set(load_transcript_data(open(fn, 'rb'), sample=sample), use_all=use_all)
        setname = os.path.basename(fn).replace('_TRANSCRIPTDATA.pickled', '').replace('Sample_', '')
        sets.append((setname, set_))
        setnames.append(setname)
    
    outfile = '+'.join(setnames)# + '.venn.txt'
    
    if not use_all:
        outfile += '.mobile_only'
        # outfile.replace('.venn.txt', '.mobile_only.venn.txt')
    
    fo = open(outfile + '.venn.txt', 'w')
    
    
    fo.write('#SETS = %i\n' % len(sets))
    fo.write('#\t'.join(['%s: %i' % (s[0], len(s[1])) for s in sets]) + '\n')
    
    processed = process_sets(sets)
    generate_diagram(outfile + '.png')
    
    for id_, set_ in processed:
        fo.write('%s\t%s\n' % (id_, len(set_)))
    fo.write('\n')
    
    for id_, set_ in processed:
        fo.write('%s\t%s\n' % (id_, len(set_)))
        for item in set_:
            fo.write('%s\n' % item)
    
    fo.close()
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

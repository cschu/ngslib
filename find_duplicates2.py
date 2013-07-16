#!/usr/bin/env python
'''
Created on Jun 25, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import re

from collections import defaultdict

class IDWithPointer(object):
    def __init__(self, id_, p=None):
        self.id_ = id_
        self.parent = None
        self.children = set()
        self.parent = p
        #if p is not None:
        #    p.children.add(self)
        if p is not None:
            self.setParent(p)
        pass
    def setParent(self, p):
        while True:
            if p is None or p.parent is None:
                break
            p = p.parent
        self.parent = p
        # print 'SELF-PARENT', self.parent.id_, self.id_, p.id_
        p.children.add(self)
        p.children = p.children.union(self.children)
        for c in self.children:
            # print '  SET PARENT', c.id_, p.id_
            c.setParent(p)
        self.children = set()
        # print 'P>', p.id_, [c.id_ for c in p.children]
        pass
        
    pass

def main(argv):

    all_objects, seq_d = {}, {}
    fn = argv[0]
    for line in open(fn):
        if line.startswith('@'):
            continue
        fields = line.strip().split()
        
        samflag = int(fields[1])
        is_mapped = samflag & 4 == 0 
                 
        if not is_mapped:
            continue
        
        p1 = sid1 = fields[0]
        p2 = sid2 = fields[2].strip('@')
        seq_d[sid1] = fields[9]
        
        allowed_match = re.search('NM:i:0', line)
        if allowed_match:
            if p2 in all_objects:
                all_objects[p1] = IDWithPointer(p1, all_objects[p2])            
            elif p1 in all_objects:
                all_objects[p2] = IDWithPointer(p2, None)
                all_objects[p1].setParent(all_objects[p2])             
            else:
                if p1 == p2:
                    all_objects[p1] = IDWithPointer(p1, None)
                else:
                    all_objects[p2] = IDWithPointer(p2, None)            
                    all_objects[p1] = IDWithPointer(p1, all_objects[p2])
        else:
            all_objects[p1] = IDWithPointer(p1, None)
        
    for k in sorted(all_objects, key=lambda x:len(all_objects[x].children), reverse=True):
        if len(all_objects[k].children) > 0 or all_objects[k].parent is None:
            print '>', all_objects[k].id_, len(all_objects[k].children)
            try:
                print seq_d[k]
            except:
                print '<SEQ MISSING>'
            for c in all_objects[k].children:
                print '\t', c.id_
                
    pass

def main2(argv):
    """
        if sid2 in maps_to:
                # if sid1 -> sid2 and sid2 -> sidX, then sid1 -> sidX
                maps_to[sid1] = find_final_target(maps_to, sid2)
            else:
                # sid1 -> sid2
                maps_to[sid1] = find_final_target(maps_to, sid2)
                for sid in maps_to:
                    if maps_to[sid] == sid1:
                        maps_to[sid] = maps_to[sid1]
    """
    
    NO_PARENT = IDWithPointer("no", None)
    
    all_objects = {}
    for p1, p2 in [('s1', 's2'), ('s3', 's1'), ('s4', 's2'), ('s2', 'sn'), ('s5', 's2'), ('s6', 's6'), ('s7', 's6')]:
        print 'Adding', p1, p2
        if p2 in all_objects:
            all_objects[p1] = IDWithPointer(p1, all_objects[p2])            
        elif p1 in all_objects:
            print '2)'
            all_objects[p2] = IDWithPointer(p2, None)
            all_objects[p1].setParent(all_objects[p2])             
        else:
            if p1 == p2:
                all_objects[p1] = IDWithPointer(p1, None)
            else:
                all_objects[p2] = IDWithPointer(p2, None)            
                all_objects[p1] = IDWithPointer(p1, all_objects[p2])
            
            
                        
        for k in sorted(all_objects):
            
            try:
                parent = all_objects[k].parent.id_
            except:
                parent = 'NO PARENT'
            print k, parent, [(c.id_, c.parent.id_) for c in all_objects[k].children]
            
        
    
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

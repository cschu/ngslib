#!/usr/bin/env python
'''
Created on Jul 12, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import time, datetime

from collections import defaultdict

def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')

            

            


class NATrie(object):
    def __init__(self):
        self.values = [] 
        self.child_A, self.child_C, self.child_G, self.child_U = None, None, None, None
        pass
    def add(self, seq, value):
        head, tail = seq[0], seq[1:]
        child = getattr(self, 'child_' + head)
        if child is None:            
            setattr(self, 'child_' + head, NATrie())
            child = getattr(self, 'child_' + head)
        if len(seq) == 1:
            child.values.append(value)
        else:
            child.add(tail, value)
        pass
    def lookup(self, seq):
        head, tail = seq[0], seq[1:]
        child = getattr(self, 'child_' + head) 
        if child is not None:
            return child.values if len(tail) == 0 else child.lookup(tail)
        return False
    def remove(self, seq):
        pass
    def seqs(self, seq=''):
        if len(self.values) > 0:
            yield seq, self.values
        for c in 'ACGU':
            child = getattr(self, 'child_' + c)
            if child is not None:
                for seq_, values in child.seqs(seq=seq + c):
                    yield seq_, values  
    pass

class DuplicateReadFinder(object):
    def __init__(self):
        self.trie = NATrie()
        self.seqd = defaultdict(list)        
        pass    
    def add_sequence(self, sid, seq):
        if 'N' in seq:
            return False
        seq = seq.replace('T', 'U')
        intree = self.trie.lookup(seq)
        if intree:            
            key = intree[0]
        else:
            key = len(self.seqd)
            self.trie.add(seq, key)            
        self.seqd[key].append(sid)
        return True          
    def write_sequences(self, out=sys.stdout):      
        all_seqs = sorted([(id_[0], len(self.seqd[id_[0]]), seq) 
                           for seq, id_ in self.trie.seqs()],
                          key=lambda x:x[1], reverse=True)
        for id_, n, seq in all_seqs:
            out.write('>%i_%i_%s\n%s\n' % (id_, n, self.seqd[id_][0], seq))            
        pass
        
        
        
        pass
    pass




def main(argv):
    trie = NATrie()    
    seqd = defaultdict(list)
    fn = argv[0]
        
    lct = 0
    sys.stderr.write('%s> Building trie...\n' % get_timestamp())
    for line in open(fn):
        if lct % 4 == 0:
            id_ = line.strip().split()[0]
        elif lct % 4 == 1:
            if not 'N' in line:                
                seq = line.strip().replace('T', 'U')
                # print id_, lct, seq                        
                intree = trie.lookup(seq)
                if intree:
                    seqd[intree[0]].append(id_)
                else:
                    trie.add(seq, len(seqd))
                    seqd[len(seqd)].append(id_)
        lct += 1
        
    sys.stderr.write('%s> Sorting sequences...\n' % get_timestamp())
    all_seqs = sorted([(id_[0], len(seqd[id_[0]]), seq) for seq, id_ in trie.seqs()],
                      key=lambda x:x[1], reverse=True)
    
    sys.stderr.write('%s> Writing sequences...\n' % get_timestamp())    
    fo = sys.stdout    
    for id_, n, seq in all_seqs:
        fo.write('%i_%i_%s\n%s\n' % (id_, n, seqd[id_][0], seq))
    # fo.close()


def main2(argv):
    
    
    trie = NATrie()
    print trie.lookup('A')
    trie.add('A', 1)
    print trie.lookup('A')
    trie.add('ACGA', 2)
    print trie.lookup('ACGA')
    trie.add('G', 3)
    trie.add('G', 4)
    print trie.lookup('G')
    print trie.lookup('AUGA')
    
    
    for s in trie.seqs():
        print s
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

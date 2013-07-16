#!/usr/bin/env python
'''
Created on Jul 4, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import datetime
import time

from collections import defaultdict


def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
 
class Trie:
    def __init__(self):
        self.root = defaultdict(Trie)
        self.value = None

    def add(self, s, value):
        """Add the string `s` to the 
        `Trie` and map it to the given value."""
        head, tail = s[0], s[1:]
        cur_node = self.root[head]
        if not tail:
            cur_node.value = value
            return  # No further recursion
        self.root[head].add(tail, value)

    def lookup(self, s, default=None):
        """Look up the value corresponding to 
        the string `s`."""
        head, tail = s[0], s[1:]
        node = self.root[head]
        if tail:
            return node.lookup(tail)
        return node.value or default

    def remove(self, s):
        """Remove the string s from the Trie. 
        Returns *True* if the string was a member.
        This leaves all prefixes of s in the Trie, which
        might cause wrong .prefix()-results.
        """
        head, tail = s[0], s[1:]
        if head not in self.root:
            return False  # Not contained
        node = self.root[head]
        if tail:
            return node.remove(tail)
        else:
            del node
            return True

    def prefix(self, s):
        """Check whether the string `s` is a prefix 
        of some member. Don't expand the trie on negatives 
        (just like the lookup-function, which also doesn't expand anything)"""
        if not s:
            # if this is left as is, it will always return None if s ends at 
            # a node that doesn't have a key ending at it
            # better: return True (otherwise, checking for 'if not trie.prefix(s)'
            # will yield a wrong result
            # return self.value
            return True
        head, tail = s[0], s[1:]
        if head not in self.root:
            return False  # Not contained
        node = self.root[head]
        return node.prefix(tail)

    def items(self):
        """Return an iterator over the items of the `Trie`."""
        for char, node in self.root.items():
            if node.value is None:
                yield node.items()
            else:
                yield node
                
    def strings(self, string=''):
        """ use dfs / preorder """
        strings = []
        if self.value is not None:
            # print 'XXX', self.value, string
            strings.append((self.value, string))
        for node in self.root:
            # print string + node, self.value
            strings.extend(self.root[node].strings(string=string + node))
        return strings
            
        
        

def main(argv):
    trie = Trie()    
    seqd = defaultdict(list)
    fn = argv[0]
        
    lct = 0
    sys.stderr.write('%s> Building trie...\n' % get_timestamp())
    for line in open(fn):
        if lct % 4 == 0:
            id_ = line.strip().split()[0]
        elif lct % 4 == 1:
            seq = line.strip()            
            intree = trie.lookup(seq)
            if intree is not None:
                seqd[intree].append(id_)
            else:
                trie.add(seq, len(seqd))
                seqd[len(seqd)].append(id_)
        lct += 1
    
    
    sys.stderr.write('%s> Sorting sequences...\n' % get_timestamp())
    all_seqs = sorted([(id_, len(seqd[id_]), seq) for id_, seq in trie.strings()],
                      key=lambda x:x[1], reverse=True)
    
    sys.stderr.write('%s> Writing sequences...\n' % get_timestamp())
    #fo = open(fn.replace('.fastq', '') + '.nodupes.fa', 'wb')
    fo = sys.stdout
    #for id_, seq in trie.strings():
    for id_, n, seq in all_seqs:
        fo.write('%i_%i\n%s\n' % (id_, n, seq))
    #fo.close()
    
    
    
    





def main_test(argv):
    
    trie = Trie()
    trie.add("GAGA", 1)
    trie.add("GUCG", 2)
    trie.add("ACGT", 4)
    trie.add("ACGTTCTGTAC", 3)
    
    
    for value, string in trie.strings():
        print value, string
    
    print trie.prefix("ACGTT")
    print trie.lookup("CUC")
    
    #items = trie.items()
    
    #for item in items:
        # print item, item.value, item.root
    #    print item
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

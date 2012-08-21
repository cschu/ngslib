#!/usr/bin/python

import os
import sys
import math

#
class CodingDictionary(dict):
    def __init__(self, open_gff, coding_flags):
        sys.stderr.write('Generating coding dictionary...\n')
        self.coding_flags = coding_flags 
        for line in open_gff.readlines():
            line = [col.strip() 
                    for col in line.split('\t')]
            gfftype = line[2]
            if gfftype in self.coding_flags.keys():
                region = xrange(int(line[3]), int(line[4]) + 1)
                src = line[0][-1] 
                self[src] = self.get(src, {})
                state = self.coding_flags[gfftype]
                for pos in region:                    
                    self[src][pos] = self[src].get(pos, 0) | state

        for src in self.keys():
            self[src] = dict(sorted(self[src].items()))
        pass
    def query_position(self, source, pos, state):
        if source in self and pos in self[source]:
            #return (self[source][pos] | state) != 0
            return (self[source][pos] & state) != 0 # QUESTION: WHICH RESULTS/METHODS DEPEND ON DATA CREATED WITH "|" ?
        return False
    pass

###
def do_something():
    return None 

###
def main(argv):
    return None

if __name__ == '__main__': main(sys.argv[1:])

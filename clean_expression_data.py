#!/usr/bin/python

'''
Created on Aug 7, 2012

@author: schudoma
'''

import sys

def main(argv):
    
    data = [line.strip() for line in open(argv[0]).readlines()]

    allow_line = True
    for line in data:
        is_agi_line = line.startswith('At')
        if is_agi_line:            
            agis = line.split()[0].split('|')
            if len(agis) > 1:
                line = line.split('\t')
                for agi in agis:
                    print '\t'.join([agi] + line[1:])
            else:
                print line
            allow_line = False
        elif allow_line:
            print line
    
    pass

if __name__ == '__main__': main(sys.argv[1:])
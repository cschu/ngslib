#!/usr/bin/env python
'''
Created on Jun 13, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import datetime
import time


def get_timestamp(file_friendly=False):
    if file_friendly:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    else:
        return datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')



def binary_search(key, lst):
    pass

def main(argv):
    
    #sys.stderr.write('# ' +  get_timestamp() + ' Reading reference list...' + '\n')
    sys.stderr.write('# ' +  get_timestamp() + ' Building reference list...' + '\n')
    sys.stderr.flush()
    reflist = set()
    last_line = ''
    for line in open(argv[0]):
        if not line.startswith('#'):
            line = line.strip().split()[1].strip('@')
            # print line
            # sys.exit()
            if line != last_line:
                last_line = line
                reflist.add(line)

    #reflist = [line.strip().strip('@').split()[1] 
    #           for line in open(argv[0])
    #           if not line.startswith('#')]
    #sys.stderr.write('# ' +  get_timestamp() + ' Removing duplicates...' + '\n')     
    #sys.stderr.flush()
    #reflist = list(set(reflist))
    # print get_timestamp(), 'Sorting...'
    # reflist.sort()
    sys.stderr.write('REFLISTLENGTH= ' + str(len(reflist)) +  '\n')
    
    
    sys.stderr.write('# ' + get_timestamp() + ' Starting search...' + '\n')
    sys.stderr.flush()
    for line in open(argv[1]):
        if not line.startswith('#'):
            #read_found = (binary_search(line.strip(), reflist) != -1)
            read_found = line.strip() in reflist
            if read_found:
                sys.stdout.write(line.strip() + '\n')
                sys.stdout.flush()

    pass

if __name__ == '__main__': main(sys.argv[1:])

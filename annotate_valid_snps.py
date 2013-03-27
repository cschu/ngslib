#!/usr/bin/env python
'''
Created on Mar 20, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys



def main(argv):
    
    valid_d = {}
    for line in open(argv[1]):
        line = line.strip().replace(' ', '\t').split('\t')
        valid_d[(line[0], line[1])] = 'valid' in line
        
    # print valid_d.keys()        
    
    for line in open(argv[0]):
        line = line.strip()
        fields = line.split('\t')
        # print line
        key = (fields[0].lstrip('Chr'), fields[10])
        if key in valid_d:
            if valid_d[key]:
                line += ';RELIABLE=YES'
            else:
                line += ';RELIABLE=NO'
            sys.stdout.write(line + '\n')
        else:
            sys.stderr.write('KEY:' + str(key) + 'NOT IN valid_d\n') 
            sys.exit(1)
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

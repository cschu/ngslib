#!/usr/bin/env python

'''
Created on Aug 21, 2012

@author: Chris
'''

import sys
import csv


class IntergenicHit(object):
    def __init__(self, fields, values):
        for field, value in zip(fields, values):
            setattr(self, field, value)            
        pass
    


def main(argv):
    
    ihits = []
    
    reader = csv.reader(open(argv[0], 'rb'), delimiter=';')
    headers = None
    for row in reader:
        if headers is None:
            headers = row
        else:
            ihit = IntergenicHit(headers, row)
            ihits.append(ihit)
    pass

if __name__ == '__main__': main(sys.argv[1:])
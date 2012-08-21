#!/usr/bin/env python

'''
Created on Aug 21, 2012

@author: Chris
'''

import sys
import csv


class IntergenicHit(object):
    def __init__(self, fields, values, casts):
        for field, value, cast in zip(fields, values, casts):
            setattr(self, field, value)            
        # print self.__dict__
        pass
    def is_valid(self, min_reads=3):
        print self.__dict__
        return self.total_reads > min_reads
    pass
    


def main(argv):
    
    obj_headers = ['contig', 'pos', 'colbase', 'pedbase', 'total_reads', 'support_ped', 'fr_ped', 'support_col', 'fr_col']
    casts = [str, int, str, str, int, int, float, int, float]
    
    reader = csv.reader(open(argv[0], 'rb'), delimiter=';')
    writer = csv.writer(open(argv[0].replace('.csv', '.filtered.csv'), 'wb'), delimiter=';', quotechar='"')
    headers = None
    for row in reader:
        if headers is None:
            headers = row
            writer.writerow(headers)
        elif row[0][0] == '#':
            writer.writerow(row)
        else:
            ihit = IntergenicHit(obj_headers, row, casts)
            if ihit.is_valid():
                # ihits.append(ihit)
                writer.writerow(row)
    
                
    pass

if __name__ == '__main__': main(sys.argv[1:])

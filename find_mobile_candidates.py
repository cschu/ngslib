#!/usr/bin/env python
'''
Created on Mar 11, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import csv


def main(argv):
    
    reader = csv.reader(open(argv[0]), delimiter=',', quotechar='"')
    writer = csv.writer(open(argv[0].replace('.csv', '_with_mobile.csv'), 'wb'), 
                        delimiter=',', quotechar='"')
    cutoff = 3
    for row in reader:
        if row[-1] == 'N/A':
            continue
        if not row[0].strip() == 'AGI':
            agi, numbers, ann = row[0], map(int, row[2:-1]), row[-1]

            mobile = []            
            if numbers[1] >= 3 or numbers[4] >= 3 or numbers[7]:
                mobile.append('ped')
            if numbers[9] >= 3 or numbers[12] >= 3 or numbers[15] >=3:
                mobile.append('col')
            
            if len(mobile) == 0:              
                continue
            else:
                row.append('/'.join(mobile))
            pass
        writer.writerow(row)
        
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

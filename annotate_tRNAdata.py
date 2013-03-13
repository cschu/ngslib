#!/usr/bin/env python
'''
Created on Mar 11, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import csv
import os

def main(argv):
    
    delimiter = '\t'
    quotechar = '"'
    
    reader = csv.reader(open(argv[0]), delimiter=delimiter, quotechar=quotechar)
    writer = csv.writer(open(argv[0].replace('.csv', '_with_%s.csv' % os.path.basename(argv[1])), 'wb'), 
                        delimiter=delimiter, quotechar=quotechar)

    annotation = {}
    for line in open(argv[1]):
        line = line.strip().split(',')
        keys = line[-2].split(';')
        val = line[-1]
        for k in keys:
            annotation[k] = val
    print annotation
    print 'AT1G71697.1' in annotation
    print 'AT1G71710.1' in annotation
    print annotation.get('AT1G71710.1', '***')
            
    for row in reader:
        if not row[0].strip() == 'AGI':
            id_ = row[17].lstrip('Parent=')
            print id_
            # if id_ == 'AT1G71710.1': 
            #    print 'XX'
            #    break 
            # print row
            # print id_
            # print annotation.get(id_, '***')
            
            row.append(annotation.get(id_, ''))
            # print row
            pass
        writer.writerow(row)        

    
    pass


def main2(argv):
    
    delimiter = ','
    quotechar = '"'
    
    reader = csv.reader(open(argv[0]), delimiter=delimiter, quotechar=quotechar)
    writer = csv.writer(open(argv[0].replace('.csv', '_with_%s.csv' % argv[1]), 'wb'), 
                        delimiter=delimiter, quotechar=quotechar)
    
    annotation = map(str.strip, open(argv[1]).readlines())
    #annotation = dict([(agi[:-2], agi) for agi in annotation])
    annotation_d = {}
    for agi in annotation:
        annotation_d[agi[:-2]] = list(set(annotation_d.get(agi[:-2], []) + [agi]))
    annotation = annotation_d
    
    for row in reader:
        if row[0].strip() == 'AGI':
            writer.writerow(row)
            continue
        # id_ = row[8].split(';')[0].lstrip('ID=')
        id_ = row[0].strip()
        
        
        if id_ in annotation:
            row.append(';'.join(annotation[id_]))
        else:
            row.append('N/A')
        writer.writerow(row)
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python

'''
Created on Sep 3, 2012

@author: schudoma
'''
import sys
import glob
import os

def read_reads(open_fn):
    last = 0
    reads = {}
    for line in open_fn:
        if line.startswith('>'):
            last = line.strip()[1:]
        else:
            reads[last] = line.strip()
    return reads

def read_bowtie(open_fn):
    reads = []
    for line in open_fn:
        line = line.strip().split('\t')
        reads.append(line[0].strip())
    return reads


def get_nomapping(all_reads, mapped_reads, output_fmt=dict):
    all_ = set(all_reads.keys())
    mapped_ = set(mapped_reads)
    nomapping = [(key, all_reads[key]) for key in list(all_ - mapped_)]
    return output_fmt(nomapping)


def main(argv):
    
    ngspath = '/home/schudoma/ngstemp/cschudoma/Sly-Cus_parasitic_RNA-Seq_Kragler'
    path_to_results = os.path.join(ngspath, 'bowtie-CTS-TTS-results')
    path_to_reads = os.path.join(ngspath, 'fasta')
    output_path = os.path.join(ngspath, 'CTS-TTS')
    
    read_files = glob.glob(os.path.join(path_to_reads, 'CTS*paired.fas'))
    for fn in read_files:
        short_fn = os.path.basename(fn.rstrip('.fas'))
        reads = read_reads(open(fn))
        print short_fn
        all_nomapping = []
        
        result_files = glob.glob(os.path.join(path_to_results, '%s_*.map' % short_fn))
        for rfn in result_files:
            results = read_bowtie(open(rfn))
            nomapping = get_nomapping(reads, results, output_fmt=list)
            all_nomapping += nomapping
            print '\t', rfn
        
        
        all_nomapping = list(dict(all_nomapping).items())
        
        print len(all_nomapping)
        
        out = open(os.path.join(output_path, short_fn + '.nomapping'), 'wb')                   
        for item in sorted(all_nomapping):
            # print item
            out.write('>%s\n%s\n' % item)
        out.close()        
        print
        # break
    
    # print read_files
                           
    
    
    
    
    
    pass



    
def main2(argv):
    all_reads = read_reads(open(argv[0]))
    mapped_reads = read_bowtie(open(argv[1]))

    out = sys.stdout
    
    for k, v in sorted(get_nomapping(all_reads, mapped_reads).items()):
        out.write('>%s\n%s\n' % (k, v))
    
    
    pass


if __name__ == '__main__': main(sys.argv[1:])
#!/usr/bin/env python
'''
Created on Jan 25, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

from get_introns_from_gff import get_attributes

def read_transcript_info(open_fn):
    transcript_info = {}
    for line in open_fn:
        line = line.strip().split(',')
        agi, is_mobile = line[0].strip('"'), line[1] == '1'
        if agi in transcript_info and transcript_info[agi] != is_mobile:
            sys.stderr.write('Conflict: %s %s vs. %s\n' % (agi, is_mobile, transcript_info[agi]))
            sys.exit(1)
        else:
            transcript_info[agi] = is_mobile
    return transcript_info

def process_intersected_gff(open_fn, transcript_info, out=sys.stdout):
    
    header = ['Contig',
              'AGI', 'tRNA_located_in', 'Start_feature', 'End_feature', 'Strand',  
              'AGI_tRNA', 'Start_tRNA', 'End_tRNA', 'Strand_tRNA', 'Length_tRNA',
              'Overlap[nt]', '%covered_tRNA', 'Transcript_found', 'Transcript_isMobile',
              'Same_strand']
    out.write(','.join(header) + '\n')    
    for line in open_fn:
        line = line.strip().split('\t')
        
        contig = line[0]
        start, end = map(int, line[3:5])
        strand = line[6]
        attributes = get_attributes(line[8])
        type_region = line[11]
        start_region, end_region = map(int, line[12:14])
        strand_region = line[15]
        attributes_region = get_attributes(line[17])
        overlap = int(line[18])
        
        trna_length = end - start + 1
        agi = attributes_region['Parent']
        if agi[:-2] in transcript_info:
            transcript_found = 'yes'
            if transcript_info[agi[:-2]]:
                transcript_is_mobile = 'yes'
            else:
                transcript_is_mobile = 'no'                
        else:
            transcript_found = 'no'
            transcript_is_mobile = 'unknown'        
            
        if strand == strand_region:
            same_strand = 'yes'
        else:
            same_strand = 'no'
        
        row = [contig, 
               agi, type_region, start_region, end_region, strand_region,
               attributes['ID'], start, end, strand, trna_length, 
               overlap, int(float(overlap)/trna_length * 100 + 0.5), 
               transcript_found, transcript_is_mobile, same_strand]
        
        out.write(','.join(map(str, row)) + '\n')
    return None
        
        
        


def main(argv):
    
    transcript_info = read_transcript_info(open(argv[0]))
    process_intersected_gff(open(argv[1]), transcript_info, out=sys.stdout)
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

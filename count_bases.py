#!/usr/bin/env python
'''
Created on Mar 14, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys


#
def count_bases(col, cutoff=-5, mult_counts=None):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 
              'low_qual': 0, 'del': 0, 'mates_disagree': 0,
              'unique_reads': 0}
    reads_at_pos = dict(zip(counts.keys(), [[] for k in counts.keys()]))  
    weighted = dict([(k, 0) for k in counts])
    bad_reads = []
    unique_reads = set([])
    valid_reads = []
    reads = {}
    count_ = 0
    
    # get all mates    
    for read in col.pileups:
        if read.alignment.flag & 0x4:
            continue
        unique_reads.add(read.alignment.qname)        
        reads[read.alignment.qname] = reads.get(read.alignment.qname, []) + [read]
    
    # check positions        
    for qname, pair in reads.items():        
        if len(pair) == 1:
            read = pair[0]            
        else:
            if pair[0].is_del == 1 and pair[1].is_del == 1:            
                bad_reads.append(('del', qname))
                weighted['del'] += 1.0 
                counts['del'] += 1
                reads_at_pos['del'].append(qname)
                continue
            elif pair[0].is_del == 1: 
                read = pair[1]
            elif pair[1].is_del == 1: 
                read = pair[0]
            else:
                if pair[0].alignment.seq[pair[0].qpos] == pair[1].alignment.seq[pair[1].qpos]: 
                    if ord(pair[0].alignment.qual[pair[0].qpos]) >= ord(pair[1].alignment.qual[pair[1].qpos]):
                        read = pair[0]
                    else:
                        read = pair[1]                        
                else:
                    # the mates do not agree at this position and have equal quality => discard read
                    weighted['mates_disagree'] += 1.0
                    counts['mates_disagree'] += 1
                    reads_at_pos['mates_disagree'].append(qname)
                    bad_reads.append(('mates_disagree', qname))
                    continue # <- important!
                pass
        base, qual = read.alignment.seq[read.qpos], ord(read.alignment.qual[read.qpos]) - 33       
        if read.is_del == 1:
            bad_reads.append(('del', qname))
            weighted['del'] += 1.0
            counts['del'] += 1
            reads_at_pos['del'].append(qname)            
        elif qual > cutoff and base != 'N':
            weighted[base] += 1.0
            counts[base] += 1
            reads_at_pos[base].append(qname)
            seq = None
            if read.indel == 0:
                seq = read.alignment.seq
            valid_reads.append((qname, base, seq, read.qpos)) 
        elif qual < cutoff:
            bad_reads.append(('low_qual', qname))
            weighted['low_qual'] += 1.0
            counts['low_qual'] += 1
            reads_at_pos['low_qual'].append(qname)
        else:
            bad_reads.append(('N', qname))
            weighted['N'] += 1.0  
            counts['N'] += 1
            reads_at_pos['N'].append(qname)
        pass
    
    counts['T'] += counts['U']
    reads_at_pos['T'].extend(reads_at_pos['U'])
    del counts['U']
    del reads_at_pos['U']
    weighted['T'] += weighted['U']
    del weighted['U']
    counts['unique_reads'] = len(unique_reads)
    reads_at_pos['unique_reads'] = list(unique_reads)

    return counts, weighted, valid_reads, reads_at_pos




#
def count_bases_old(col, cutoff=-5, mult_counts=None):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0, 
              'low_qual': 0, 'del': 0, 'mates_disagree': 0}
    weighted = dict([(k, 0) for k in counts])
    bad_reads = []
    unique_reads = set([])
    valid_reads = []
    reads = {}
    count_ = 0
    
    # get all mates    
    for read in col.pileups:
        if read.alignment.flag & 0x4:
            continue
        unique_reads.add(read.alignment.qname)        
        reads[read.alignment.qname] = reads.get(read.alignment.qname, []) + [read]
    
    # check positions        
    for qname, pair in reads.items():
        # read = pair[0]
        
        # base1, base2 = pair[0].alignment.seq[pair[0].qpos], None
        # qual1, qual2 = ord(pair[0].alignment.qual[pair[0].qpos]) - 33, None
        # isdel1, isdel2 = pair[0].is_del == 1, None        
        
        if len(pair) == 1:
            # base, qual, isdel = base1, qual1, isdel1
            read = pair[0]            
        else:
            # base2 = pair[1].alignment.seq[pair[1].qpos]
            # qual2 = ord(pair[1].alignment.qual[pair[1].qpos]) - 33
            # isdel2 = pair[1].is_del == 1
            
            # if isdel1 and isdel2:
            if pair[0].is_del == 1 and pair[1].is_del == 1:            
                bad_reads.append(('del', qname))
                weighted['del'] += 1.0 
                counts['del'] += 1
                continue
            elif pair[0].is_del == 1: #isdel1:
                # base, qual, isdel = base2, qual2, False
                read = pair[1]
            elif pair[1].is_del == 1: #isdel2:
                # base, qual, isdel = base1, qual1, False
                read = pair[0]
            else:
                # if base1 == base2:
                if pair[0].alignment.seq[pair[0].qpos] == pair[1].alignment.seq[pair[1].qpos]: 
                    # mates agree => check base  
                    # base, qual, isdel = base1, max(qual1, qual2), False
                    if ord(pair[0].alignment.qual[pair[0].qpos]) >= ord(pair[1].alignment.qual[pair[1].qpos]):
                        read = pair[0]
                    else:
                        read = pair[1]                        
                else:
                    # the mates do not agree at this position and have equal quality => discard read
                    weighted['mates_disagree'] += 1.0
                    counts['mates_disagree'] += 1
                    bad_reads.append(('mates_disagree', qname))
                    continue # <- important!
                pass
        # if isdel:
        base, qual = read.alignment.seq[read.qpos], ord(read.alignment.qual[read.qpos]) - 33       
        if read.is_del == 1:
            bad_reads.append(('del', qname))
            weighted['del'] += 1.0
            counts['del'] += 1            
        elif qual > cutoff and base != 'N':
            # print 'x', base, qual
            weighted[base] += 1.0
            counts[base] += 1
            seq = None
            if read.indel == 0:
                seq = read.alignment.seq
            valid_reads.append((qname, base, seq, read.qpos)) #, read.alignment.pos, read.alignment.aend))
        elif qual < cutoff:
            bad_reads.append(('low_qual', qname))
            weighted['low_qual'] += 1.0
            counts['low_qual'] += 1
        else:
            bad_reads.append(('N', qname))
            weighted['N'] += 1.0  
            counts['N'] += 1
        pass
    
    counts['T'] += counts['U']
    del counts['U']
    weighted['T'] += weighted['U']
    del weighted['U']
    counts['unique_reads'] = len(unique_reads)

    return counts, weighted, valid_reads




def main(argv):
    pass

if __name__ == '__main__': main(sys.argv[1:])

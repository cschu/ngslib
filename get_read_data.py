#!/usr/bin/env python
'''
Created on Mar 1, 2013

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys

"""
 mit den Col/Ped SNPs und den experimentell gefundenen Reads. 
 Wenn Du mir fuer ein Transkript die entsprechenden SNPs und die Reads 
 (mit Sequenz, Anfangsposition, Laenge, Natur der polymorphen Nukleotide) 
"""

def main(argv):
    
    lines = []
    for line in open(argv[0]):
        line = line.strip().split('\t')
        lines.append(line)
        # break 
    lines = sorted(lines, key=lambda x:(x[0], int(x[1])))
    ['HWI-ST377:127:D0PHGACXX:7:1102:16868:187523', '0', 'Chr3', '20794121', '40', '16M77N92M259N64M', '*', '0', '0', 'TATGTATCCAGTTTGTCTTCTAGACTGCCGTGGCATCTAGATGGAAACTTATTAAAAAAGTTGTCTGCTATGTTTTGCATCGTGCATTCCAATCTACCACCCATAACGCTGTGATAGTCTCCATACTGATCAAGGTACTCAATTCGCTTCTTTGTATTGAGAACCATATTTG', 'CCCFFFFFHHHHHJJIIJJJJJIJJJJIJJIIIJJJJJJJJJJJJJJJJJJJJJJIJJJJFHIIIJJIJHHHHHHFFFFFFDDDDDEEDEEEDDDDDDDFHHHHIHJJJJJJJJJJJJJIIHJIIJJJJJJJJJIHEJJJJJJJJJJJJJJJJJJJJJJHHHGHFFFFFCCC', 'MD:Z:172', 'NH:i:1', 'HI:i:1', 'NM:i:0', 'SM:i:40', 'XQ:i:40', 'X2:i:0', 'XS:A:-', 'RG:Z:col-wt-root.R1.extended.fq.gsnap.unpaired_uniq']

    read_d = {}
    for line in open(argv[1]):
        line = line.strip().split('\t')
        read_d[line[1]] = line[2]    
    
    
    header = ['ReadID', 'Supported_Accession', 'RNA_insert_length', 
              'Start_Mate1', 'Seq_Mate1', 'Length_Mate1', 
              'Start_Mate2', 'Seq_Mate2', 'Length_Mate2']
    print '\t'.join(header)
    
    current = None
    for line in lines:
        read_data = [line[0], line[1], line[3], line[5], line[8], line[9], len(line[9])]
        #if 'paired_no' in line[-1]:
        if int(read_data[1]) & 2 == 2:            
            if current is None:
                current = read_data
            else:
                if current[0] != read_data[0]:
                    print 'PROBLEM'
                    print current
                    print line
                    sys.exit()
                    
                #rdata = [current[0], '%i/%i' % (current[1], read_data[1]), current[2], 
                #         '%s+%s' % (current[3], read_data[3]), ]
                # print 'PE1\t' + '\t'.join(map(str, current))
                # print 'PE2\t' + '\t'.join(map(str, read_data))
                
                positions = map(int, [current[2], read_data[2]])
                if min(positions) == positions[0]:
                    # 'current' is first mate
                    seq1, len1, seq2, len2 = current[5], current[6], read_data[5], read_data[6]
                else:  
                    seq1, len1, seq2, len2 = read_data[5], read_data[6], current[5], current[6] 
                
                rdata = [current[0], read_d[current[0]], abs(int(current[4])), 
                         min(positions), seq1, len1, max(positions), seq2, len2]
                print '\t'.join(map(str, rdata))
                
                current = None
        else:
            # del read_data[1]
            rdata = [read_data[0], read_d[read_data[0]], read_data[6], read_data[2], 
                     read_data[5], read_data[6], '', '', '']
            
            
            print '\t'.join(map(str, rdata))
                
    
    
    # paired = sorted([line for line in lines if 'paired_no' in line[-1]])
    # unpaired = [line for line in lines if not 'paired_no' in line[-1]]
    
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

#!/usr/bin/env python
'''
Created on Aug 24, 2012

@author: schudoma
'''
import re
import sys

# ported from https://bitbucket.org/galaxy/galaxy-central/src/tip/tools/samtools/pileup_parser.pl

def count_bases_in_pileup(bases, qualities, refbase, qual_cutoff=-5):
    # counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'bad': 0}
    counts = {}
    
    bases = re.sub('\^.', '', bases)
    bases = re.sub('$', '', bases)
    
    """
     if ($read_bases =~ m/[\$\^\+-]/) {
        $read_bases =~ s/\^.//g; #removing the start of the read segement mark
        $read_bases =~ s/\$//g; #removing end of the read segment mark
        while ($read_bases =~ m/[\+-]{1}(\d+)/g) {
            my $indel_len = $1;
            $read_bases =~ s/[\+-]{1}$indel_len.{$indel_len}//; # remove indel info from read base field
        }
    }
    """
    while True:
        mobj = re.search('[\+-][0-9]+', bases)
        if not mobj:
            break
        # print mobj.group()
        indel_len = int(mobj.group()[1:])
        # print '[\+-]%i[ACGTNacgtn]{%i}' % (indel_len, indel_len)
        bases = re.sub('[\+-]%i[ACGTNacgtn]{%i}' % (indel_len, indel_len), 
                       '', bases)
        # print bases, len(bases)
        # break
        pass

    bases_above_qual_cutoff = 0    
    for base, qual in zip(list(bases), list(qualities)):
        if (ord(qual) - 33 >= qual_cutoff) and base != '*':            
            if base in 'acgtACGT':
                counts[base] = counts.get(base, 0) + 1
                bases_above_qual_cutoff += 1  
            elif base == '.' or base == ',':
                counts[refbase] = counts.get(refbase, 0) + 1
                bases_above_qual_cutoff += 1  
            elif base == '>' or base == '<':
                counts['skip'] = counts.get('skip', 0) + 1
            else:
                counts['weird'] = counts.get('weird', []) + [base]
        else:
            counts['bad'] = counts.get('bad', 0) + 1
    
    return counts

def parse_pileupfile(lines, offset=4):
    # for line in open_file:
    p = 0
    while p < len(lines):
        
        
        pos, pileup = lines[p:p+offset][0], lines[p:p+offset][offset - 1]
        # Chr1    3218123 A       11      tttTTTTTTT,     HFJJJJJHJJF
        # print pos
        # print pileup
        pileup = pileup.split()
        
        if len(pileup) >= 5:
            read_bases, read_qualities = pileup[4:]
            counts = count_bases_in_pileup(read_bases, read_qualities, pos.split(':')[-1].strip()[0])
            if len(counts) > 0:
                print '=======>', pos, pileup
                print pos, pileup, counts
        p += offset
        # print p
    pass
        


def main(argv):
    
    parse_pileupfile(open(argv[0]).readlines())
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])


"""
'ref skips' in mpileup-output: ignore the read 

Chr1:2738-2740:CT
Chr1    2738    T    2    >>    :G
Chr1    2739    A    2    >>    :G
Chr1    2740    C    2    >>    :G

samtools view col-N-root.all.sorted.bam Chr1:2738-2740
HWI-ST377:146:D0YLAACXX:7:1311:9412:39906       163     Chr1    2581    39      8M2850N63M    =5507    129     AAGAAAAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCT 8?;B?DB3:2A+A?GHCFA;AADFGCGE@)1:?;H@HE@B?8B=8=)32@;DCE===7A=);9?>>3(;>3        MD:Z:71 NH:i:1  HI:i:1NM:i:0   SM:i:39 XQ:i:40 X2:i:0  XS:A:+  RG:Z:col-N-root.R1.paired_no.fq.gsnap.concordant_uniq
HWI-ST377:146:D0YLAACXX:7:2302:18885:20969      163     Chr1    2581    40      8M2850N92M    =5551    220     AAGAAAAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTG   BCCFDDFFG?CFHIIIFHIIIJIJJJIGIIJGGIIJGJFFHGHJJIGIJJG>DHI@CHHHFFEFCEDEEEDDACCCBACDBDDBDDB2<BCDCCC3<?@C   MD:Z:100        NH:i:1  HI:i:1  NM:i:0  SM:i:40 XQ:i:40 X2:i:0XS:A:+   RG:Z:col-N-root.R1.paired_no.fq.gsnap.concordant_uniq


"""



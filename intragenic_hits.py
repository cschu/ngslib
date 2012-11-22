#!/usr/bin/env python
'''
Created on Nov 20, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import sys
import os

from collections import defaultdict

"""
header = ['AGI', 'Pos_SNP', '#Reads', '#Reads_Col', '#Reads_Ped', 
          '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']
"""
header = ['#Reads', '#Reads_Col', '#Reads_Ped', 
          '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']

class AthGene(object):
    def __init__(self, agi, snp=None):
        self.agi = agi
        self.snps = {snp.pos: snp}
        pass
    pass

class AthSNP(object):
    def __init__(self, pos):
        self.pos = pos
        self.readcounts = defaultdict(dict)        
        pass
    def add_readcount(self, sampleID, counts, keys=header):
        try:
            self.readcounts[sampleID] = dict(zip(keys, counts))
        except:
            sys.stderr.write('Could not process readcounts:\nkeys:%s\ncounts%s' % (keys, counts))
            sys.exit(1)
            
        # print '>>', sampleID, self.readcounts[sampleID]#, condition
        pass
    
    def __repr__(self):         
        return ',\n'.join(map(str, self.readcounts))
    def get_counts_string(self, counts, sampleID, delimiter=','):
        return delimiter.join([str(self.readcounts[sampleID][count]) 
                               for count in counts])
    
    pass

"""
header = ['AGI', 'Pos_SNP', '#Reads', '#Reads_Col', '#Reads_Ped', 
          '#Reads_N', '#Reads_lowqual', '#Reads_del', '#Reads_other']
"""

def write_perGene_summary(genes, sampleIDs, count_mask, header, delimiter=',', out=sys.stdout):
    
    out.write('%s\n' % delimiter.join(header))
    for agi, snps in sorted(genes.items()):        
        snp_counts = dict(zip(sampleIDs, 
                              [dict(zip(count_mask + ['coverage'], [0, 0, 0, 0])) for sid in sampleIDs]))        
        row = [agi, ''] 
        for pos, counts in sorted(snps.items()):            
            # row = [agi] #, '', len(snps)] # AGI, N_SNPS (needed from external source), N_coveredSNPS
            for sid in sorted(sampleIDs):
                if sid in counts.readcounts:
                    for cm in count_mask:
                        snp_counts[sid][cm] += counts.readcounts[sid][cm]
                    if counts.readcounts[sid]['#Reads_Col'] > 0 or counts.readcounts[sid]['#Reads_Ped'] > 0:
                        snp_counts[sid]['coverage'] += 1 
        for sid in sampleIDs:
            row.extend(map(str, [snp_counts[sid][cm] for cm in count_mask]))
            row.append(str(snp_counts[sid]['coverage']))
              
        out.write('%s\n' % delimiter.join(row))        
    pass


def write_perSNP_summary(genes, sampleIDs, count_mask, header, delimiter=',', out=sys.stdout):
    
    out.write('%s\n' % delimiter.join(header))    
    for agi, snps in sorted(genes.items()):
        for pos, counts in sorted(snps.items()):
            row = [agi, str(pos)]
            for sid in sorted(sampleIDs):
                if sid in counts.readcounts:
                    row.extend(counts.get_counts_string(count_mask,sid).split(','))   
                else:
                    row.extend(['0' for item in count_mask])
            out.write('%s\n' % delimiter.join(row))
    pass
    


def main(argv):
    
    args = argv    
    genes = {}    
    sampleIDs = set()
    out = sys.stdout
    delimiter = ','
    
    for fn in args:
        open_fn = open(fn)
        fn = os.path.basename(fn).strip()
        sampleID = fn[:fn.find('.')]
        sampleIDs.add(sampleID)
        
        for line in open_fn:            
            line = line.strip().split(',')
            if line[0] == 'AGI':
                continue
            snp_pos = int(line[1])
            if line[0] not in genes:
                genes[line[0]] = {snp_pos: AthSNP(snp_pos)}
            gene = genes[line[0]]                        
            if snp_pos not in gene:
                gene[snp_pos] = AthSNP(snp_pos)            
            gene[snp_pos].add_readcount(sampleID, map(int, line[2:]))

    sampleIDs = sorted(list(sampleIDs))               
    count_mask = ['#Reads_Col', '#Reads_Ped', '#Reads']
    out_header = [('%s:Col,%s:Ped,%s:Total_Reads' % (sid, sid, sid)).split(',')
                  for sid in sampleIDs]
    out_header = ['AGI', 'SNP'] + [item for sublist in out_header for item in sublist]
    
    write_perSNP_summary(genes, sampleIDs, count_mask, out_header, delimiter=delimiter, 
                         out=open('intragenic_snp_summary.csv', 'wb'))  
    
    out_header = [('%s:Col,%s:Ped,%s:Total_Reads,%s:SNPs_covered' % (sid, sid, sid, sid)).split(',')
                  for sid in sampleIDs]
    out_header = ['AGI', '#SNPs'] + [item for sublist in out_header for item in sublist]
    #print out_header
    write_perGene_summary(genes, sampleIDs, count_mask, out_header, delimiter=delimiter, 
                          out=open('intragenic_gene_summary.csv', 'wb'))
    
    
       
    pass

if __name__ == '__main__': main(sys.argv[1:])

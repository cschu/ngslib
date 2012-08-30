#!/usr/bin/env python

'''
Created on Aug 15, 2012

@author: Chris
'''
import sys

try:
    import pysam
except:
   pass


from snptools import SNPDict
import gff_helpers
import find_covered_genes as FIND_GENES


MIN_NREADS = 3


# region = (gffline[0], int(gffline[3]) - 1, 
#                  int(gffline[4]) - 1, gffline[6], gene_id)

# self[(snpline.seqid, snpline.pos)] = snpline Chr%

#
def find_region(pos, regions):
    """
    uses binary search adaptation from 
    http://stackoverflow.com/a/212413/902449
    """
    low, high = 0, len(regions)
    while low < high:
        mid = (low + high) // 2
        region = regions[mid]
        print pos, region, low, high,
        if pos[0] < region[0]:
            print 'RULE 1'
            # position is located on a contig 
            # with lesser id than current region
            high = mid
        elif pos[0] > region[0]:
            print 'RULE 2'
            # position is located on a contig 
            # with higher id than current region
            low = mid + 1
        else:
            # position is located on same contig as region ...
            if pos[1] < region[1]:
                print 'RULE 3'
                # ... but upstream of region
                high = mid
            elif pos[1] > region[2]:
                print 'RULE 4'
                # ... but downstream of region
                low = mid + 1
            else:
                print 'RULE 5'
                # and is also located within region
                return region
    print 'RULE 6'
    return None

###
def remove_intragenic_snps(snp_d, intragenic_regions):
    for pos in sorted(snp_d.keys()):
        #pos = ('Chr1', 5936)
        region = find_region(pos, intragenic_regions)        
        if not region is None:
            # print 'REMOVING', pos, snp_d[pos], region
            del snp_d[pos]
        #sys.exit()
        
    
    return snp_d

def main(argv):
    intragenic_regions = gff_helpers.read_intragenic_regions(open('tair10_genes.gff'))
    intragenic_regions = sorted(intragenic_regions)
    
    # x = find_region(('Chr3', 16092855), intragenic_regions)
    # print x
    
    x = find_region(('Chr1', 26179018), intragenic_regions)
    print x
    
    #snp_d = SNPDict(open('ped-0-snps_no-indels.txt'))
    #snp_d = remove_intragenic_snps(snp_d, intragenic_regions)
    
    #print ('Chr3', 16092855) in snp_d
    #print ('Chr1', 26179018) in snp_d
    
    #for k, v in sorted(snp_d.items()):
    #    print k, v
    #pass

###
def main2(argv):

    samfile = pysam.Samfile(argv[0], 'rb')
    fo = sys.stdout
    # fo = open('%s.covered_genes.csv' % argv[0].rstrip('.bam'), 'w')    
    # fo2 = open('%s.covered_genes_with_reads.csv' % argv[0].rstrip('.bam'), 'w')
    fo2=sys.stdout

    # fo.write('%s\n' % ','.join(COL_HEADERS))
    
    """
    # tair10_genes.gff
    intragenic_regions = gff_helpers.read_intragenic_regions(open(argv[1]))
    intragenic_regions = sorted(intragenic_regions)
    # print intragenic_regions[:10]
    
    # ped-0-snps_no-indels.txt
    snp_d = SNPDict(open(argv[2]))  
    # print snp_d.items()[:10]
    snp_d = remove_intragenic_snps(snp_d, intragenic_regions)
    """
    snps = gff_helpers.read_snp_from_gff(open(argv[1]))
    
    count_snps = 0
    print ';'.join(['contig', 'position', 'refbase_(Col)', 'mutation_(Ped)', 'total_reads', 
                    '#support_Ped', 'fr_Ped', '#support_Col', 'fr_Col'])
    
    # for snp_id, snpline in sorted(snp_d.items()):
    for snp_id in snps:        
        basecount = FIND_GENES.count_bases(samfile, snp_id[0], snp_id[1])
    
        refbase = basecount.get(snp_id[3], 0.0)
        snpbase = basecount.get(snp_id[4], 0.0)
        
        total_reads = sum(basecount.values()) - basecount['bad']
        
        if total_reads > 0:
            
            # line = str(snpline).split('\t')[1:5]
            
            # snp = (gffline[0], int(gffline[3] + 1), int(gffline[4]), comments['refbase'], comments['mutation'])
            
            line = [snp_id[0], snp_id[1], snp_id[3], snp_id[4]]            
            
            line.extend([total_reads, 
                         snpbase, float(snpbase)/total_reads,
                         refbase, float(refbase)/total_reads])
            print ';'.join(map(str, line))
                    
            # print snpline, 'x', total_reads, snpbase, float(snpbase)/total_reads,
            # print refbase, float(refbase)/total_reads 
            count_snps += 1
    print '# Total SNPs: %i Covered: %i (%.3f)' % (len(snp_d), count_snps, float(count_snps)/len(snp_d)) 
    
    # tair10_genes.gff
    # process_gff(open(argv[3]), polymorphs, snp_d, samfile, fo, fo2, min_reads=MIN_NREADS)
    # fo2.close()
    
    
    fo.close()
    samfile.close()
    return None



if __name__ == '__main__': main(sys.argv[1:])

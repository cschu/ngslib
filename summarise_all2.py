#!/usr/bin/env python
'''
Created on Nov 30, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import os
import sys
import csv
import glob
import pickle

SNP_CUTOFF = 3

def process_hits(open_fn):
    # print open_fn.name
    reader = csv.reader(open_fn, delimiter=',', quotechar='"')
    total = 0
    
    hits = {}
    for row in reader:
        
        region = row[0].strip()
        # print region
        if region in ('AGI', 'Contig'):
            header = row
            continue
        else:
            row_d = dict(zip(header[2:], map(float, row[2:-2])))
        
        if region not in hits:            
            hits[region] = {}
        if row[1] not in hits[region]:
            hits[region][row[1]] = {'col': 0, 'ped': 0, 'other': 0}
        for k1, k2 in [('#Reads_Col', 'col'), ('#Reads_Ped', 'ped'), ('#Reads_other', 'other')]:
            hits[region][row[1]][k2] = hits[region][row[1]].get(k2, 0) + row_d[k1]
            total += row_d[k1]
            
    return hits, total

def normalise_RPM(C, N):
    return 10**6 * C / float(N)

def find_mobileRNAs(hits_d, sampleIDs, normalise={}, cutoff=3, out=sys.stdout):
    
    mobile_counts = dict(zip(sampleIDs, [set([]) for sid in sampleIDs]))
    for region, region_d in sorted(hits_d.items()):        
        for snp, snp_d in sorted(region_d.items()):   
            for sampleID in sampleIDs:
                if sampleID in snp_d:
                    col, ped = snp_d[sampleID]['col'], snp_d[sampleID]['ped']
                    if col >= cutoff and 'ped' in sampleID:
                        mobile_counts[sampleID].add(region)                         
                    elif ped >= cutoff and 'col' in sampleID:
                        mobile_counts[sampleID].add(region)
    for sid in sorted(sampleIDs):
        out.write('%s\t%i\t%f\n' % (sid, len(mobile_counts[sid]), 
                                    normalise_RPM(len(mobile_counts[sid]), normalise.get(sid, 10**6))))
    pass

def find_mobileRNAs_nodupes(rpg, normalise={}, cutoff=3, out=sys.stdout):
    
    mobile_counts = dict(zip(rpg.keys(), [0 for sid in rpg.keys()]))
    for sid, agis in rpg.items():
        for agi, counts in agis.items():
            col, ped = counts['col'], counts['ped']
            if col >= cutoff and 'ped' in sid:
                mobile_counts[sid] += 1
            elif ped >= cutoff and 'col' in sid:
                mobile_counts[sid] += 1
        pass
    for sid in rpg:
        out.write('%s\t%i\t%f\n' % (sid, mobile_counts[sid], 
                                    normalise_RPM(mobile_counts[sid], normalise.get(sid, 10**6))))
    pass


def filter_regions(hits_d, snp_cutoff=SNP_CUTOFF):
    for region, region_d in hits_d.items():
        for snp, snp_d in region_d.items():
            for sampleID, sample_d in snp_d.items():
                if sample_d['col'] < snp_cutoff and sample_d['ped'] < snp_cutoff:
                    del hits_d[region][snp][sampleID]
            if len(hits_d[region][snp]) == 0:
                del hits_d[region][snp]
        if len(hits_d[region]) == 0:
            del hits_d[region]
    return hits_d
                    
                         

def write_perSNP_summary(hits_d, sampleIDs, col0, normalise=None, out=sys.stdout):
    header = [col0, 'SNP'] + ['%s:col,%s:ped,%s:hits' % (sid, sid, sid)
                              for sid in sampleIDs]
    out.write('%s\n' % ','.join(header))
    for region, region_d in sorted(hits_d.items()):        
        for snp, snp_d in sorted(region_d.items()):            
            row = [region, snp]
            for sampleID in sampleIDs:                 
                if sampleID in snp_d:
                    values = [snp_d[sampleID]['col'],
                              snp_d[sampleID]['ped'],
                              sum(snp_d[sampleID].values())]
                    if normalise is not None:
                        values = map(lambda x:normalise_RPM(x, normalise[sampleID]), values)
                    row.extend(values)
                else:
                    row.extend([0, 0, 0])
            out.write('%s\n' % ','.join(map(str, row)))
    pass

def write_perGene_summary(hits_d, sampleIDs, normalise=None, out=sys.stdout):
    header = ['AGI', '#SNPs'] + ['%s:col,%s:ped,%s:hits,%s:snp-covered' % (sid, sid, sid, sid)
                                 for sid in sampleIDs] 
    out.write('%s\n' % ','.join(header)) 
    transcript_count = dict(zip(sampleIDs, [0 for sid in sampleIDs]))
    for region, region_d in sorted(hits_d.items()):
        region_counts = dict(zip(sampleIDs, [{'col': 0, 'ped': 0, 'other': 0, 'cov': 0} 
                                             for sid in sampleIDs]))  
        # print region_counts      
        for snp, snp_d in sorted(region_d.items()):
            for sampleID in sampleIDs:                 
                if sampleID in snp_d:
                    region_counts[sampleID]['cov'] += 1
                    for k, v in snp_d[sampleID].items():                        
                        region_counts[sampleID][k] += v
        row = [region, '']
        for sampleID in sampleIDs:
            values = [region_counts[sampleID]['col'],
                      region_counts[sampleID]['ped'],
                      sum(region_counts[sampleID].values()) - region_counts[sampleID]['cov']]
            if sum(region_counts[sampleID].values()) - region_counts[sampleID]['cov'] > 0:
                transcript_count[sampleID] += 1
            if normalise is not None:
                values = map(lambda x:normalise_RPM(x, normalise[sampleID]), values)
            row.extend(values)
            row.append(region_counts[sampleID]['cov'])
        out.write('%s\n' % ','.join(map(str, row)))
    return transcript_count

# nodupes means that reads that hit two or more SNPs are only counted once            
def write_perGeneSummary_nodupes(rpg, sampleIDs, normalise=None, out=sys.stdout, cutoff=3):
    header = ['AGI', '#SNPs'] + ['%s:col,%s:ped,%s:hits' % (sid, sid, sid)
                                 for sid in sampleIDs] 
    out.write('%s\n' % ','.join(header))
    # this works even if stupid eclipse shows an error!
    
    regions = sorted(set([item for sublist in [v.keys() for v in rpg.values()] for item in sublist]))
    transcript_count = dict(zip(sampleIDs, [0 for sid in sampleIDs]))
    for region in regions:
        # if region != 'AT1G33370': continue
        row = [region, '']
        allow_region = False
        for sampleID in sampleIDs:
            values = [0, 0, 0]
            if region in rpg[sampleID]:
                # print region, sampleID, 'XXX'
                
                total = rpg[sampleID][region]['col'] + rpg[sampleID][region]['ped']
                # 2013-01-21: why is the sum of col and ped reads used as cutoff criterion?
                # if total >= cutoff:
                if rpg[sampleID][region]['col'] >= cutoff or rpg[sampleID][region]['ped'] >= cutoff:
                    allow_region = True
                    values = [rpg[sampleID][region]['col'],
                              rpg[sampleID][region]['ped'],
                              total]
                    transcript_count[sampleID] += 1
                # print values, max(values), cutoff                                
                # values.append(sum(values))
                # if max(values[:2]) >= cutoff:                    
                if normalise is not None:
                    values = map(lambda x:normalise_RPM(x, normalise[sampleID]), values)
            row.extend(values)
        if allow_region:
            out.write('%s\n' % ','.join(map(str, row)))
    return transcript_count
            
                          
                          


def main(argv):
    
    args = argv
    if 'intra' in args[0]:
        col0 = 'AGI'
        mode = 'intra'
    else:
        col0 = 'Contig'
        mode = 'inter'
    
    rpg = {}
    total_reads = {}    
    for fn in glob.glob('*.%s.rpg.dat' % mode):
        obj = pickle.load(open(fn, 'rb'))
        # for k, v in
        sid = fn[:fn.find('.')]
        print fn, sid        
        if sid not in rpg:
            print 'new sid:', sid
            total_reads[sid] = sum([counts['col'] + counts['ped'] 
                                    for counts in obj.values()])
            rpg[sid] = obj
        else:
            for agi, counts in obj.items():
                total_reads[sid] += counts['col'] + counts['ped']
                if agi not in rpg[sid]:
                    rpg[sid][agi] = {}
                for k in counts:
                    rpg[sid][agi][k] = rpg[sid][agi].get(k, 0) + counts[k]                
        pass
            
    for k, c in total_reads.items():
        print '%s: %i reads\n' % (k, c)
    for k, c in rpg.items():
        print '%s: %i genes\n' % (k, len(c))
    # return None
    
    """
    try:
        prefix = args[0][:args[0].find('.') + 1]
    except:
        prefix = '002.'
    """
    prefix = '003.'
    
    all_hits = {}
    sampleIDs = set([])
    
    # total_reads = {}
    for fn in args:
        sampleID = os.path.basename(fn)
        sampleID = sampleID.lstrip(prefix)
        sampleID = sampleID[:sampleID.find('.')]
        sampleIDs.add(sampleID)        
        
        hits_d, total = process_hits(open(fn))
        ## TODO
        # 0. total reads needs to be obtained from rpg -- maybe total reads needs to come from raw count?        
        # 1. write_perGene_summary needs to use rpg-data for the correct read/gene count!
        # total_reads[sampleID] = total_reads.get(sampleID, 0) + total        
        for region, region_d in hits_d.items():
            all_hits[region] = all_hits.get(region, {})
            for snp, snp_d in region_d.items():            
                all_hits[region][snp] = all_hits[region].get(snp, {})
                all_hits[region][snp][sampleID] = all_hits[region][snp].get(sampleID, dict(zip(snp_d.keys(), 
                                                                                               [0 for i in xrange(len(snp_d))])))
                for count, value in snp_d.items():
                    all_hits[region][snp][sampleID][count] += value
    
    # throw away any low-populated snps and genes
    all_hits = filter_regions(all_hits)
    sampleIDs = sorted(list(sampleIDs))
    
    
    
    write_perSNP_summary(all_hits, sampleIDs, col0, out=open('%sgenic_snp_summary.csv' % mode, 'wb'))
    write_perSNP_summary(all_hits, sampleIDs, col0, 
                         normalise=total_reads, out=open('%sgenic_snp_summary_normalised.csv' % mode, 'wb'))
    if mode == 'intra':
        tcount = write_perGene_summary(all_hits, sampleIDs, out=open('intragenic_gene_summary.csv', 'wb'))
        write_perGene_summary(all_hits, sampleIDs, 
                              normalise=total_reads, out=open('intragenic_gene_summary_normalised.csv', 'wb'))
        out = open('transcript_count.txt', 'wb')
        for sid in sorted(tcount.keys()):
            out.write('%s\t%i\n' % (sid, tcount[sid]))
        out.close()
        find_mobileRNAs(all_hits, sampleIDs, normalise=total_reads, out=open('mobileRNAs.txt', 'wb'))
        
        tcount = write_perGeneSummary_nodupes(rpg, sampleIDs, out=open('intragenic_gene_summary_nodupes.csv', 'wb'))
        write_perGeneSummary_nodupes(rpg, sampleIDs,
                                     normalise=total_reads,
                                     out=open('intragenic_gene_summary_nodupes_normalised.csv', 'wb'))
        out = open('transcript_count_nodupes.txt', 'wb')
        for sid in sorted(tcount.keys()):
            out.write('%s\t%i\n' % (sid, tcount[sid]))
        out.close()
        find_mobileRNAs_nodupes(rpg, normalise=total_reads, out=open('mobileRNAs_nodupes.txt', 'wb'))
        
        
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

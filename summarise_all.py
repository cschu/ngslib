#!/usr/bin/env python
'''
Created on Nov 30, 2012

@author: Christian Schudoma (schudoma@mpimp-golm.mpg.de, cschu@darkjade.net)
'''
import os
import sys
import csv

SNP_CUTOFF = 3

def process_hits(open_fn):
    reader = csv.reader(open_fn, delimiter=',', quotechar='"')
    total = 0
    
    hits = {}
    for row in reader:
        region = row[0].strip()
        if region in ('AGI', 'Contig'):
            header = row
            continue
        else:
            row_d = dict(zip(header[2:], map(float, row[2:])))
        
        # if row_d['#Reads_Col'] >= SNP_CUTOFF or row_d['#Reads_Ped'] >= SNP_CUTOFF:
        if region not in hits:
            # hits[region] = hits.get(region, {row[1]: {'col': 0, 'ped': 0, 'other': 0}})
            hits[region] = {row[1]: {'col': 0, 'ped': 0, 'other': 0}}
        for k1, k2 in [('#Reads_Col', 'col'), ('#Reads_Ped', 'ped'), ('#Reads_other', 'other')]:
            # hits[region][row[1]][k2] += row_d[k1]
            hits[region][row[1]][k2] = hits[region][row[1]].get(k2, 0) + row_d[k1]
            total += row_d[k1]
            # hits[region][row[1]] = {'col': row_d['#Reads_Col'], 
            #                        'ped': row_d['#Reads_Ped'],
            #                        'other': row_d['#Reads_other']}
            # total += row_d['#Reads_other']
            # total += sum(hits[region][row[1]].values())
            
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
                    if col > cutoff and 'ped' in sampleID:
                        mobile_counts[sampleID].add(region)                         
                    elif ped > cutoff and 'col' in sampleID:
                        mobile_counts[sampleID].add(region)
    for sid in sorted(sampleIDs):
        out.write('%s\t%i\t%f\n' % (sid, len(mobile_counts[sid]), 
                                   normalise_RPM(len(mobile_counts[sid]), normalise.get(sid, 10**6))))
    pass
                  
def filter_regions(hits_d, snp_cutoff=SNP_CUTOFF):
    for region, region_d in hits_d.items():
        for snp, snp_d in region_d.items():
            for sampleID, sample_d in snp_d:
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
            



def main(argv):
    
    args = argv
    if 'intra' in args[0]:
        col0 = 'AGI'
        mode = 'intra'
    else:
        col0 = 'Contig'
        mode = 'inter'
    
    try:
        prefix = args[0][:args[0].find('.') + 1]
    except:
        prefix = '000.'
    
    all_hits = {}
    sampleIDs = set([])
    
    total_reads = {}
    for fn in args:
        sampleID = os.path.basename(fn)
        sampleID = sampleID.lstrip(prefix)
        sampleID = sampleID[:sampleID.find('.')]
        sampleIDs.add(sampleID)        
        
        hits_d, total = process_hits(open(fn))
        total_reads[sampleID] = total_reads.get(sampleID, 0) + total        
        for region, region_d in hits_d.items():
            all_hits[region] = all_hits.get(region, {})
            for snp, snp_d in region_d.items():            
                all_hits[region][snp] = all_hits[region].get(snp, {})
                all_hits[region][snp][sampleID] = all_hits[region][snp].get(sampleID, dict(zip(snp_d.keys(), 
                                                                                               [0 for i in xrange(len(snp_d))])))
                for count, value in snp_d.items():
                    all_hits[region][snp][sampleID][count] += value
    # print all_hits
    # print
    # print
    all_hits = filter_regions(all_hits)
    
    write_perSNP_summary(all_hits, sorted(list(sampleIDs)), col0, out=open('%sgenic_snp_summary.csv' % mode, 'wb'))
    write_perSNP_summary(all_hits, sorted(list(sampleIDs)), col0, 
                         normalise=total_reads, out=open('%sgenic_snp_summary_normalised.csv' % mode, 'wb'))
    if mode == 'intra':
        tcount = write_perGene_summary(all_hits, sorted(list(sampleIDs)), out=open('intragenic_gene_summary.csv', 'wb'))
        write_perGene_summary(all_hits, sorted(list(sampleIDs)), 
                              normalise=total_reads, out=open('intragenic_gene_summary_normalised.csv', 'wb'))
        out = open('transcript_count.txt', 'wb')
        for sid in sorted(tcount.keys()):
            out.write('%s\t%i\n' % (sid, tcount[sid]))
        out.close()
        
        find_mobileRNAs(all_hits, sorted(list(sampleIDs)), normalise=total_reads, out=open('mobileRNAs.txt', 'wb'))
    
    
    pass

if __name__ == '__main__': main(sys.argv[1:])

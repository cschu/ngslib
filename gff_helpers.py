'''
Created on Jul 11, 2012

@author: schudoma
'''

#
def parse_gff_comments(string):
    return dict([x.split('=')
                 for x in string.strip().split(';')])

#
def read_polymorphs(open_gff):
    polymorphs = {}
    for gffline in open_gff:
        gffline = gffline.strip().split('\t')
        comments = parse_gff_comments(gffline[8])
        snp = (gffline[0], int(gffline[3]) - 1, 
               int(gffline[4]) - 1, gffline[6])
        gene_id = comments['ID']
        polymorphs[gene_id] = polymorphs.get(gene_id, []) + [snp]
    return polymorphs

#
def read_snp_from_gff(open_gff):
    snps = []
    for gffline in open_gff:
        gffline = gffline.strip().split('\t')
        comments = parse_gff_comments(gffline[8])
        snp = (gffline[0], int(gffline[3] + 1), int(gffline[4]), comments['refbase'], comments['mutation'])
        snps.append(snp)
    return snps


def read_intragenic_regions(open_gff):
    i_regions = []
    for gffline in open_gff:
        gffline = gffline.strip().split('\t')
        comments = parse_gff_comments(gffline[8])
        gene_id = comments['ID']
        region = (gffline[0], int(gffline[3]) - 1, 
                  int(gffline[4]) - 1, gffline[6], gene_id)
        i_regions.append(region)
    return i_regions

#
def get_geneid(gffline):
    comments_d = parse_gff_comments(gffline[8])     
    return gffline[:1] + gffline[3:5] + [comments_d.get('ID', 'N/A'),
                                         comments_d.get('Note', 'N/A')]

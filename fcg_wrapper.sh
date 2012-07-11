#!/bin/bash

NGSPATH=/home/schudoma/projects/ngs
POLYMORPH_GFF=tair10_ped-0_polymorphic_genes.gff
SNP_DICT=ped-0-snps_no-indels.txt
GENES_GFF=tair10_genes.gff

FCG_PY=/home/schudoma/workspace/ngslib/find_covered_genes.py

python $FCG_PY $1 $NGSPATH/$POLYMORPH_GFF $NGSPATH/$SNP_DICT $NGSPATH/$GENES_GFF
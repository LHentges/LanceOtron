#!/usr/bin/env python
# coding: utf-8

import numpy as np
import csv
import re

motif_pfm = 'MA0139.1.pfm'
genome_fasta_file = 'hg38.fa'
out_file = '../motif_enrichment/ChIP-seq_CTCF_spleen/CTCF_exact-match_genome.bed'


nucleotide_freq_list = []
with open(motif_pfm) as f:
    next(f)
    row_reader = csv.reader(f, delimiter=' ')
    for row in row_reader:
        float_list = [float(x) for x in row]
        nucleotide_freq_list.append(float_list)

nucleotide_freq_array = np.array(nucleotide_freq_list).T

nucleotide_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
motif_nuc = ''
for row in nucleotide_freq_array:
    most_common_nuc = row.argmax()
    if row[most_common_nuc]/row.sum()>=(.75):
        motif_nuc+=nucleotide_dict[most_common_nuc]
    else:
        motif_nuc+='.'

motif_nuc_rev = ''
nucleotide_dict_rev = {'A':'T', 'C':'G', 'G':'C', 'T':'A', '.':'.'}
for nuc in motif_nuc:
    motif_nuc_rev = nucleotide_dict_rev[nuc]+motif_nuc_rev    

motif_bed_list = []
with open(genome_fasta_file) as f:
    for row in f:
        if row.startswith('>'):
            chrom = row[1:-1]
            chrom_row_counter = 0
        for m in re.finditer(motif_nuc, row.upper()):
            motif_bed_list.append([chrom, (chrom_row_counter*50)+m.start(), (chrom_row_counter*50)+m.end()])
        for m in re.finditer(motif_nuc_rev, row.upper()):
            motif_bed_list.append([chrom, (chrom_row_counter*50)+m.start(), (chrom_row_counter*50)+m.end()])
        chrom_row_counter+=1

with open(out_file, 'w', newline='') as f:
    bed_writer = csv.writer(f, delimiter='\t')
    bed_writer.writerows(motif_bed_list)
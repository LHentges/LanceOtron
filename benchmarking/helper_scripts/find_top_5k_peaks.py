#!/usr/bin/env python
# coding: utf-8

import pybedtools
import csv
import os

sample_list = ['ChIP-seq_H3K27ac_HAP1', 'ChIP-seq_H3K4me3_MG63']
peakcaller_type_list = ['0-LanceOtron', '1-MACS2', '2-LanceOtron-with-input', '3-MACS2-with-input']
out_folder = '../TSS_annotations/'

def bed_file_to_list(bed_file, header=False):
    bed_list = []
    with open(bed_file, 'r') as f:
        bed_reader = csv.reader(f, delimiter='\t')
        if header:
            next(bed_reader)
        for row in bed_reader:
            row_list = []
            for i, column in enumerate(row):
                if (i==1) or (i==2):
                    row_list.append(int(column))
                else:
                    row_list.append(column)
            bed_list.append(row_list)
    return bed_list


# sorting on column 4 for LanceOtron, overall peak score, and column 9 for MACS2, q-value

for sample in sample_list:
    for peakcaller_type in peakcaller_type_list:
        bed_file = '../labelled_datasets/'+sample+'/'+peakcaller_type+'_'+sample
        if peakcaller_type.startswith('0') or peakcaller_type.startswith('2'):
            bed_file+='.bed'
            bed_list = bed_file_to_list(bed_file)
            bed_list.sort(key=lambda x: x[3], reverse=True)
        elif peakcaller_type.startswith('1') or peakcaller_type.startswith('3'):
            bed_file+='.narrowPeak'
            bed_list = bed_file_to_list(bed_file)
            bed_list.sort(key=lambda x: x[8], reverse=True)
        bed_list_centered = []
        top_peaks_count = 0
        for bed_entry in bed_list:
            if bed_entry[0]!='chrM':
                region_size = bed_entry[2]-bed_entry[1]
                mid = int((region_size)/2)+bed_entry[1]
                bed_list_centered.append([bed_entry[0], mid-500, mid+500, bed_entry[4]])
                top_peaks_count+=1
            if top_peaks_count==5000:
                break
        bed_list_centered_pybedtools = pybedtools.BedTool(bed_list_centered)
        out_file = peakcaller_type+'-top-5k-peaks_'+sample+'.bed'
        bed_list_centered_pybedtools.sort().saveas(out_folder+sample+'/'+out_file)


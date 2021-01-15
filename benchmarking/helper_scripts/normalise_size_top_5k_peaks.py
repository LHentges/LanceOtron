#!/usr/bin/env python
# coding: utf-8

import pybedtools
import csv
import os

peak_calls_folder = '../peak_calls/'

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

peak_call_list = []
with open('peak_calls_for_normalise_size_top_5k_peaks.txt') as f:
    for peak_call in f:
        peak_call_list.append(peak_call.strip())

# sorting on fourth column for LanceOtron, overall peak score, and ninth column for MACS2, q-value

for peak_call in peak_call_list:
    if peak_call.startswith('0') or peak_call.startswith('2'):
        sort_column = 3
    elif peak_call.startswith('1') or peak_call.startswith('3'):
        sort_column = 8
    bed_list = bed_file_to_list(peak_calls_folder+peak_call)
    bed_list.sort(key=lambda x: x[sort_column], reverse=True)
    bed_list_centered = []
    top_peaks_count = 0
    for bed_entry in bed_list:
        if bed_entry[0]!='chrM' and not bed_entry[0].startswith('chrUn'):
            region_size = bed_entry[2]-bed_entry[1]
            mid = int((region_size)/2)+bed_entry[1]
            bed_list_centered.append([bed_entry[0], mid-500, mid+500, bed_entry[sort_column]])
            top_peaks_count+=1
        if top_peaks_count==5000:
            break
    bed_list_centered_pybedtools = pybedtools.BedTool(bed_list_centered)
    peak_call_name, ext = peak_call.split('.')
    bed_list_centered_pybedtools.sort().saveas(peak_calls_folder+peak_call_name+'_top-5k-peaks-normalised_to_1kb.bed')


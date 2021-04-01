#!/usr/bin/env python
# coding: utf-8

import pybedtools
import os
import csv

out_folder = '../results/'
peak_calls_folder = '../peak_calls/'
peak_call_pairs_file = 'peak_call_pairs_for_find_exclusive_peaks.txt'

peak_call_pairs_list = []
with open (peak_call_pairs_file) as f:
    csv.reader = csv.reader(f)
    for peak_call_pair in csv.reader:
        peak_call_pairs_list.append(peak_call_pair)
        
for peak_call_1, peak_call_2 in peak_call_pairs_list:
    peak_call_1_name_list = peak_call_1.split('.')
    peak_call_1_detail_list = peak_call_1_name_list[0].split('_')
    peak_call_2_name_list = peak_call_2.split('.')
    peak_call_2_detail_list = peak_call_2_name_list[0].split('_')
    exclusive_peak_call_1 = peak_call_1_name_list[0]+'_'+peak_call_2_detail_list[0]+'_exclusive-peaks.'+peak_call_1_name_list[1]
    exclusive_peak_call_2 = peak_call_2_name_list[0]+'_'+peak_call_1_detail_list[0]+'_exclusive-peaks.'+peak_call_2_name_list[1]
    
    peak_call_1_pybedtools = pybedtools.BedTool(peak_calls_folder+peak_call_1)
    peak_call_2_pybedtools = pybedtools.BedTool(peak_calls_folder+peak_call_2)
    peak_call_1_intersecting_peaks = peak_call_1_pybedtools.intersect(peak_call_2_pybedtools, u=True)
    peak_call_2_intersecting_peaks = peak_call_2_pybedtools.intersect(peak_call_1_pybedtools, u=True)
    peak_call_1_exclusive_peaks = peak_call_1_pybedtools.intersect(peak_call_2_pybedtools, v=True, wa=True).saveas(peak_calls_folder+exclusive_peak_call_1)
    peak_call_2_exclusive_peaks = peak_call_2_pybedtools.intersect(peak_call_1_pybedtools, v=True, wa=True).saveas(peak_calls_folder+exclusive_peak_call_2)
    
    with open(out_folder+'exclusive-peak-calls_'+peak_call_1_detail_list[0]+'_'+peak_call_2_detail_list[0]+'.txt', 'w') as f:
        f.write(peak_call_1)
        f.write('\n--------------------------------------------------\n')
        f.write('- total peaks called: {}\n'.format(len(peak_call_1_pybedtools)))
        f.write('- intersecting peaks: {}\n'.format(len(peak_call_1_intersecting_peaks)))
        f.write('- exclusive peaks called: {}\n'.format(len(peak_call_1_exclusive_peaks)))
        f.write('\n\n')
        f.write(peak_call_2)
        f.write('\n--------------------------------------------------\n')
        f.write('- total peaks called: {}\n'.format(len(peak_call_2_pybedtools)))
        f.write('- intersecting peaks: {}\n'.format(len(peak_call_2_intersecting_peaks)))
        f.write('- exclusive peaks called: {}\n'.format(len(peak_call_2_exclusive_peaks)))
        

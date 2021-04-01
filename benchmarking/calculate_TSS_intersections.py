#!/usr/bin/env python
# coding: utf-8

import pybedtools
import os

sample_folder = './TSS_annotations/'
tss_annotation_bed_file = 'refTSS_v3.1_human_coordinate_merged.hg38.bed'
peak_calls_folder = './peak_calls/'
peak_calls_file = 'peak_calls_for_calculate_TSS-histone-ChIP-seq_intersections.txt'
out_folder = './results/'

sample_list = []
for sample in os.listdir(sample_folder):
    if os.path.isdir(sample_folder+sample):
        sample_list.append(sample)
        
sample_peak_calls_dict = {}
for sample in sample_list:
    peak_call_list = []
    with open(sample_folder+sample+'/'+peak_calls_file) as f:
        for peak_call in f:
            peak_call_list.append(peak_calls_folder+peak_call.strip())
    sample_peak_calls_dict[sample] = peak_call_list

tss_annotation_pybedtools = pybedtools.BedTool(sample_folder+tss_annotation_bed_file)
for sample in sample_list:
    with open(out_folder+sample+'_TSSs_overlapping_peaks.txt', 'w') as f:
        f.write(sample+'\n\n')
    for peak_call in sample_peak_calls_dict[sample]:
        peak_call_pybedtools = pybedtools.BedTool(peak_call)
        peaks_intersecting_tsss = peak_call_pybedtools.intersect(tss_annotation_pybedtools, u=True)
        peak_caller = peak_call.split('/')[-1].split('_')[0]
        with open(out_folder+sample+'_TSSs_overlapping_peaks.txt', 'a') as f:
            f.write(peak_caller)
            f.write('\n--------------------------------------------------\n')
            f.write('total peaks called: {}\n'.format(len(peak_call_pybedtools)))
            f.write('peaks overlapping TSSs: {}\n'.format(len(peaks_intersecting_tsss)))
            f.write('% peaks intersecting TSSs: {}\n'.format(((len(peaks_intersecting_tsss)/len(peak_call_pybedtools))*100)))
            f.write('\n\n')

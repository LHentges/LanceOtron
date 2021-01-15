#!/usr/bin/env python
# coding: utf-8

import pybedtools
import os

sample_folder = './TSS_annotations/'
out_folder = './results/'
tss_annotation_bed_file = 'refTSS_v3.1_human_coordinate_merged.hg38.bed'

sample_list = []
for sample in os.listdir(sample_folder):
    if sample!=tss_annotation_bed_file:
        sample_list.append(sample)

tss_annotation_pybedtools = pybedtools.BedTool(sample_folder+tss_annotation_bed_file)
for sample in sample_list:
    with open(out_folder+sample+'_TSSs_overlapping_top_5k_peaks.txt', 'w') as f:
        f.write(sample+'\n\n')
    for peak_call in sorted(os.listdir(sample_folder+sample)):
        peak_call_pybedtools = pybedtools.BedTool(sample_folder+sample+'/'+peak_call)
        tsss_intersecting_peaks = peak_call_pybedtools.intersect(tss_annotation_pybedtools)
        peak_caller = peak_call.split('/')[-1].split('_')[0]
        with open(out_folder+sample+'_TSSs_overlapping_top_5k_peaks.txt', 'a') as f:
            f.write(peak_caller)
            f.write('\n--------------------------------------------------\n')
            f.write('TSSs overlapping top 5k peaks: {}\n'.format(len(tsss_intersecting_peaks)))
            f.write('\n\n')
            

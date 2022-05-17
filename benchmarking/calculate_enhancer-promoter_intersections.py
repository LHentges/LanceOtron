#!/usr/bin/env python
# coding: utf-8

import pybedtools
import os
import numpy as np
import matplotlib.pyplot as plt

sample_folder = './enhancer-promoter_annotations/'
out_folder = './results/'
peak_calls_folder = './peak_calls/'
peak_calls_file = 'peak_calls_for_calculate_enhancer-promoter_intersections.txt'
enhancer_annotation_bed_file = 'GenoSTAN_enhancers.bed'
promoter_annotation_bed_file = 'GenoSTAN_promoters.bed'
use_theme_colors = True

colors_dict = {'20-LanceOtron-with-input':['darkgrey', '#C24D51'], '30-MACS2-with-input':['darkgrey', '#8172B2']}

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

enhancer_pybedtools = pybedtools.BedTool(sample_folder+enhancer_annotation_bed_file)
promoter_pybedtools = pybedtools.BedTool(sample_folder+promoter_annotation_bed_file)
for sample in sample_list:
    with open(out_folder+sample+'_peaks_intersecting_enhancers-promoters.txt', 'w') as f:
        f.write(sample+'\n\n')
    for peak_call in sorted(sample_peak_calls_dict[sample]):
        peak_call_pybedtools = pybedtools.BedTool(peak_call)
        peaks_intersecting_enhancers = peak_call_pybedtools.intersect(enhancer_pybedtools, u=True)
        peaks_not_intersecting_enhancers = peak_call_pybedtools.intersect(enhancer_pybedtools, v=True)
        peaks_intersecting_promoters = peak_call_pybedtools.intersect(promoter_pybedtools, u=True)
        peaks_not_intersecting_promoters = peak_call_pybedtools.intersect(promoter_pybedtools, v=True)
        
        not_enhancers_intersecting_promoters = peaks_not_intersecting_enhancers.intersect(promoter_pybedtools, u=True)
        #not_promoters_intersecting_enhancers = peaks_not_intersecting_promoters.intersect(enhancer_pybedtools, u=True)
        peaks_intersecting_enhancers_promoters = len(not_enhancers_intersecting_promoters)+len(peaks_intersecting_enhancers)
        peak_caller = peak_call.split('/')[-1].split('_')[0]
        with open(out_folder+sample+'_peaks_intersecting_enhancers-promoters.txt', 'a') as f:
            f.write(peak_caller)
            f.write('\n--------------------------------------------------\n')
            f.write('total peaks: {}\n'.format(len(peak_call_pybedtools)))
            f.write('\n')
            f.write('peaks intersecting enhancers: {}\n'.format(len(peaks_intersecting_enhancers)))
            f.write('peaks not intersecting enhancers: {}\n'.format(len(peaks_not_intersecting_enhancers)))
            f.write('peaks intersecting promoters: {}\n'.format(len(peaks_intersecting_promoters)))
            f.write('peaks not intersecting promoters: {}\n'.format(len(peaks_not_intersecting_promoters)))
            f.write('\n')
            f.write('peaks intersecting enhancers or promoters: {}\n'.format(len(not_enhancers_intersecting_promoters)+len(peaks_intersecting_enhancers)))
            f.write('% peaks intersecting enhancers or promoters: {}\n'.format((((peaks_intersecting_enhancers_promoters)/len(peak_call_pybedtools))*100)))
            f.write('\n\n')
        fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect='equal'))
        data = [len(peak_call_pybedtools)-peaks_intersecting_enhancers_promoters, peaks_intersecting_enhancers_promoters]
        if use_theme_colors:
            ax.pie(data, wedgeprops=dict(width=0.5), startangle=90, colors=colors_dict[peak_caller])
        else:
            ax.pie(data, wedgeprops=dict(width=0.5), startangle=90)
        ax.set_title(peak_caller)
        plt.savefig(out_folder+sample+'_'+peak_caller+'_enhancer-promoter_pie.png', bbox_inches='tight')


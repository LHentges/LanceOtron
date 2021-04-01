#!/usr/bin/env python
# coding: utf-8

import pybedtools
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sample_folder = './motif_enrichment/'
out_folder = './results/'
peak_calls_folder = './peak_calls/'
peak_calls_file = 'peak_calls_for_count_peaks_intersecting_motif.txt'
use_theme_colors = True

plt.style.use('seaborn')
colors = [["#4C72B1", "#55A868","#C24D51", '#8172B2', '#CCB974'], ['darkgrey']]

sample_list = []
for sample in os.listdir(sample_folder):
    sample_list.append(sample)

sample_peak_calls_dict = {}
sample_motif_dict = {}
for sample in sample_list:
    peak_call_list = []
    with open(sample_folder+sample+'/'+peak_calls_file) as f:
        for peak_call in f:
            peak_call_list.append(peak_calls_folder+peak_call.strip())
    sample_peak_calls_dict[sample] = peak_call_list
    for file in os.listdir(sample_folder+sample):
        if file.endswith('.bed'):
            sample_motif_dict[sample] = sample_folder+sample+'/'+file
            
for sample in sample_list:
    motif_pybedtools = pybedtools.BedTool(sample_motif_dict[sample])
    sample_results_dict = {}
    with open(out_folder+sample+'_peaks_intersecting_motifs.txt', 'w') as f:
        f.write(sample+'\n\n')
    for peak_call in sorted(sample_peak_calls_dict[sample]):
        peak_call_pybedtools = pybedtools.BedTool(peak_call)
        peaks_intersecting_motifs = peak_call_pybedtools.intersect(motif_pybedtools, u=True)
        peak_caller = peak_call.split('/')[-1].split('_')[0]
        total_peaks_count = len(peak_call_pybedtools)
        peaks_intersecting_motifs_count = len(peaks_intersecting_motifs)
        peaks_not_intersecting_motifs_count = total_peaks_count-peaks_intersecting_motifs_count
        sample_results_dict[peak_caller] = {'motif':peaks_intersecting_motifs_count, 'no motif':peaks_not_intersecting_motifs_count}
        sample_results_df = pd.DataFrame(sample_results_dict)
        with open(out_folder+sample+'_peaks_intersecting_motifs.txt', 'a') as f:
            f.write(peak_caller)
            f.write('\n--------------------------------------------------\n')
            f.write('total peaks called: {}\n'.format(total_peaks_count))
            f.write('peaks intersecting motif: {}\n'.format(peaks_intersecting_motifs_count))
            f.write('% peaks intersecting motif: {}\n'.format((peaks_intersecting_motifs_count/total_peaks_count)*100))
            f.write('\n\n')
    if use_theme_colors:
        sample_results_df.T.plot.bar(color=colors, secondary_y='no motif', legend=None)
    else:
        sample_results_df.T.plot.bar(secondary_y='no motif')
    plt.title(sample, color='black')
    plt.savefig(out_folder+sample+'_motifs_in_peak_calls.png', bbox_inches='tight')
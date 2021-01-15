#!/usr/bin/env python
# coding: utf-8

import pybedtools
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sample_folder = './labelled_datasets/'
out_folder = './results/'
plt.style.use('seaborn')

sample_list = []
for sample in os.listdir(sample_folder):
    sample_list.append(sample)

sample_file_dict = {}
for sample in sample_list:
    peak_call_file_list = []
    for file in os.listdir(sample_folder+sample):
        if file.startswith('peaks'):
            labelled_peaks_bed_file = sample_folder+sample+'/'+file
        elif file.startswith('noise'):
            labelled_peaks_noise_file = sample_folder+sample+'/'+file
        else:
            peak_call_file_list.append(sample_folder+sample+'/'+file)
    sample_file_dict[sample] = {'labelled_peaks_bed_file':labelled_peaks_bed_file, 'labelled_noise_bed_file':labelled_peaks_noise_file, 'peak_call_file_list':peak_call_file_list}

for sample in sorted(sample_file_dict):
    labelled_peaks_pybedtools = pybedtools.BedTool(sample_file_dict[sample]['labelled_peaks_bed_file'])
    labelled_noise_pybedtools = pybedtools.BedTool(sample_file_dict[sample]['labelled_noise_bed_file'])
    model_list = []
    precision_list = []
    recall_list = []
    specificity_list = []
    results_dict = {}
    with open(out_folder+sample+'_model_performance_benchmarks.txt', 'w') as f:
        f.write(sample+'\n\n')
    for peak_call in sorted(sample_file_dict[sample]['peak_call_file_list']):
        peak_call_pybedtools = pybedtools.BedTool(peak_call)
        true_positive_regions = labelled_peaks_pybedtools.intersect(peak_call_pybedtools, u=True, wa=True)
        false_negative_regions = labelled_peaks_pybedtools.intersect(peak_call_pybedtools, v=True, wa=True)
        true_negative_regions = labelled_noise_pybedtools.intersect(peak_call_pybedtools, v=True, wa=True)
        false_positive_regions = labelled_noise_pybedtools.intersect(peak_call_pybedtools, u=True, wa=True)
        precision = len(true_positive_regions)/(len(true_positive_regions)+len(false_positive_regions))
        recall = len(true_positive_regions)/(len(true_positive_regions)+len(false_negative_regions))
        specificity = len(true_negative_regions)/(len(true_negative_regions)+len(false_positive_regions))
        f1 = len(true_positive_regions)/(len(true_positive_regions)+(0.5*(len(false_positive_regions)+len(false_negative_regions))))
        peak_caller = peak_call.split('/')[-1].split('_')[0]
        results_dict[peak_caller] = {'precision':precision, 'recall/sensitivity':recall, 'specificity':specificity, 'F1':f1}
        with open(out_folder+sample+'_model_performance_benchmarks.txt', 'a') as f:
            f.write(peak_caller)
            f.write('\n--------------------------------------------------\n')
            f.write('- true positives: {}\n'.format(len(true_positive_regions)))
            f.write('- false negatives: {}\n'.format(len(false_negative_regions)))
            f.write('- true negatives: {}\n'.format(len(true_negative_regions)))
            f.write('- false positives: {}\n'.format(len(false_positive_regions)))
            f.write('\n')
            f.write('- precision: {}\n'.format(precision))
            f.write('- recall/sensitivity: {}\n'.format(recall))
            f.write('- specificity: {}\n'.format(specificity))
            f.write('- F1 score: {}\n'.format(f1))
            f.write('\n\n')
    results_df = pd.DataFrame(data=results_dict)
    results_df.plot.bar()
    plt.title(sample, color='black')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.savefig(out_folder+sample+'.png', bbox_inches='tight')

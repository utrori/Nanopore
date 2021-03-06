# coding: utf-8
import math
from visualize_nanopore_read import plot_read_structure
from ont_fast5_api.fast5_interface import get_fast5_file
import os
import search_rDNA_reads
from rDNA_structure_of_long_reads import easy_flag
import subprocess
from Bio.Seq import Seq
import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import utilities


def find_cpg(read):
    cpg_sites = []
    read = read.upper()
    for n in range(len(read)-1):
        if read[n] == 'C':
            if read[n+1] == 'G':
                cpg_sites.append(n)
    return cpg_sites


rDNA_read_ids = []
id2direction = {}
with open('reads_visualized.txt') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_id = item[0].split('@')[1].strip()
        rDNA_read_ids.append(read_id)
        summary = 0
        read_vis = item[1]
        for item in read_vis:
            if item =='￫':
                summary += 1
            elif item =='￩':
                summary -= 1
        if summary > 0:
            id2direction[read_id] = '+'
        elif summary < 0:
            id2direction[read_id] = '-'

def guppy_fast5_extraction(filename):
    data_set = []
    with get_fast5_file(filename, mode='r') as f:
        for read_id in f.get_read_ids():
            read = f.get_read(read_id)
            latest_basecall = read.get_latest_analysis('Basecall_1D')
            fastq = read.get_analysis_dataset(latest_basecall, 
            'BaseCalled_template/Fastq')
            mod_base_table = read.get_analysis_dataset(latest_basecall, 
            'BaseCalled_template/ModBaseProbs')
            data_set.append((read_id, fastq, mod_base_table))
    return data_set


def make_methylation_summary(read, mod_base_table):
    cpgs = find_cpg(read)
    mod_scores = []
    for cpg in cpgs:
        scores = mod_base_table[cpg]
        mod_scores.append((cpg, scores[3]/255))
    return mod_scores

fast5_names = os.listdir('guppy_with_mod/workspace/')
fast5_files = ['guppy_with_mod/workspace/' + i for i in fast5_names]

split_length = 200
for fast5 in fast5_files:
    dataset = guppy_fast5_extraction(fast5)
    for data in dataset:
        read_id = data[0]
        split_fastq = data[1].split()
        read = split_fastq[6]
        print(len(read))
        quality = split_fastq[8]
        mod_base_table = data[2]
        mod_scores = make_methylation_summary(read, data[2])
        utilities.split_mapping_and_sam_analysis(split_length, read_id, read, quality, 'rDNA_index/humRibosomal.fa')
        lc = plot_read_structure('test', split_length, 0)
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        x = []
        y = []
        for score in mod_scores:
            x.append(score[0])
            y.append(score[1] * 10000)
        for n, item in enumerate(mod_base_table):
            if item[1] > 120:
                ax.bar(n, item[1] * 60, width=300, color='red', zorder=3)
        ax.bar(x, y, width=100, color = 'mediumblue', zorder=0)
        ax.add_collection(lc)
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
        ax.autoscale()
        ax.set_ylim([-12000, 46000])
        plt.savefig('guppy_methylation_figs2/' + read_id + '.png', dpi=300)
        plt.close()

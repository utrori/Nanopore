# -*- coding: utf-8 -*-

import sys
import utilities
import itertools
import os
import shutil
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import numpy as np
import glob


class Fastqread(object):
    def _is_this_healthy_rDNA(self):
        """Check if the read is completely inside rDNA or at the end.
        """
        if self.length < 20000:
            return 0
        mapping_state = []
        for item in self.sam_summary:
            if item[1] != '0':
                mapping_state.append(1)
            else:
                mapping_state.append(0)
        threshold = 0.8
        if mapping_state == [] or max(mapping_state) == 0:
            return 0
        if sum(mapping_state)/len(mapping_state) > threshold:
            return 1
        else:
            for i in range(1, len(mapping_state) - 50):
                if sum(mapping_state[i:])/len(mapping_state[i:]) > threshold or \
                   sum(mapping_state[:-i])/len(mapping_state[:-i]) > threshold:
                    healthy = 2
        return 0
 
    def __init__(self, fastq_tuple, ref):
        self.seq = fastq_tuple[1].strip()
        self.quality = fastq_tuple[3].strip()
        self.read_id = fastq_tuple[0].split()[0][1:]
        self.length = len(self.seq)
        self.split_length = 300
        self.sam_summary = np.array(utilities.split_mapping_and_sam_analysis(self.split_length, 'temp', self.seq, self.quality, ref))
        self.health = self._is_this_healthy_rDNA()

    def search_invertion(self):
        directions = self.sam_summary[:,1]
        plus = len([1 for i in directions if i == '+']) + 0.1
        minus = len([1 for i in directions if i == '-']) + 0.1
        if 0.1 < plus / minus < 10:
            self.plot_structure(ref, 'temp_figs')

    def plot_structure(self, ref, savedir):
        lc = utilities.plot_read_structure(self.read_id, self.split_length, self.sam_summary)
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        ax.add_collection(lc)
        average_qual = np.mean([ord(i) - 33 for i in self.quality])
        low_quality_bases = []
        for i in range(0, len(self.quality)-20):
            if np.mean([ord(j) - 33 for j in self.quality[i:i+20]]) + 4 < average_qual:
                low_quality_bases.append(i)
        ax.bar(low_quality_bases, [10000 for i in range(len(low_quality_bases))])
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
        ax.autoscale()
        ax.set_ylim([-12000, 46000])
        if savedir[-1] == '/':
            savedir = savedir[:-1]
        plt.savefig(savedir + '/' + self.read_id + '.png', dpi=300)
        plt.close()


def find_rDNA_reads(fastq_path, output_name):
    ret_str = ''
    with open(fastq_path) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            r = Fastqread(items, 'rDNA_index/humRibosomal.fa')
            if r.health != 0:
                ret_str += ''.join(items)
    with open(output_name, 'w') as fw:
        fw.write(ret_str)


if __name__ == '__main__':
    find_rDNA_reads('/mnt/data2/hpgp/HG02717_1_Guppy_4.0.11_prom.fastq', '/mnt/data2/hpgp/HG02717_1_rDNA.fastq')
    find_rDNA_reads('/mnt/data2/hpgp/HG02717_1_Guppy_4.0.11_prom.fastq', '/mnt/data2/hpgp/HG02717_1_rDNA.fastq')
    quit()
    shutil.rmtree('temp_figs')
    os.mkdir('temp_figs')
    ref = 'rDNA_index/humRibosomal.fa'
    with open('inverted_reads.fastq') as f:
        for n, item in enumerate(itertools.zip_longest(*[iter(f)]*4)):
            r = Fastqread(item, ref)
            r.plot_structure(ref, 'temp_figs')

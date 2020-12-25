# -*- coding: utf-8 -*-

import random
import sys
#import fast5class
import utilities
import itertools
import os
import shutil
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import numpy as np
import glob


split_length = 300


class Reference(object):
    def __init__(self, ref_filename):
        self.path = ref_filename
        seq = ''
        with open(ref_filename) as f:
            self.header = f.readline().split()[0]
            for line in f:
                seq += line.strip()
        self.seq = seq


class Fastqread(object):
    def _is_this_healthy_rDNA(self):
        """Check if the read is completely inside rDNA or at the end.

        Args:
            None
        Returns:
            healthy (int): 0 if it is not healthy, 1 if it is completely inside, 
            2 if it is a read at the end of rDNA repeat.
        """
        #length cutoff
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
        self.sam_summary = np.array(utilities.split_mapping_and_sam_analysis(split_length, 'temp', self.seq, self.quality, ref.path))
        self.health = self._is_this_healthy_rDNA()

    def search_invertion(self):
        """Check if the read is inverted based on just the number of direction in split mapped reads.

        Args:
            None
        Returns:
            None (plots the read if it is inverted)
        """
        directions = self.sam_summary[:,1]
        plus = len([1 for i in directions if i == '+']) + 0.1
        minus = len([1 for i in directions if i == '-']) + 0.1
        if 0.1 < plus / minus < 10:
            self.plot_structure('temp_figs')

    def plot_structure(self, savedir):
        """Check if the read is inverted based on just the number of direction in split mapped reads.

        Args:
            savedir (str): the destination of plot
        Returns:
            None (plots the read)
        """
        lc = utilities.plot_read_structure(self.read_id, split_length, self.sam_summary)
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

    def _find_coordinate(self, coord):
        """Find exact rDNA coordinate in the given nanopore reads.

        Args:
            read (str): read sequence
            coord (int): target coordinate of rDNA

        Returns:
            result: list of [mapped coordinate, direction of mapping]
        """
        result = []
        temp_fastq_length = 500
        with open('temp_index/temp_index.fasta', 'w') as fw:
            fw.write('>{}\n{}'.format(self.header, self.seq))
        with open('temp_index/coordinate_rDNA.fastq', 'w') as fw:
            fw.write('>temp\n{}\n+\n{}\n'.format(ref.seq[coord-1:coord+temp_fastq_length-1], 'J' * temp_fastq_length))
        # with -a option, multiple hits are more clearly shown
        utilities.bwa_mapping('temp_index/temp_index.fasta', 'temp_index/coordinate_rDNA.fastq', 'temp_index/temp_sam4coord.sam')
        with open('temp_index/temp_sam4coord.sam') as samf:
            map_result = samf.readlines()[2:]
        for mapping in map_result:
            row = mapping.strip().split()
            AS = int(mapping.strip().split('AS:i:')[1].split()[0])
            flag = int(row[1])
            if easy_flag(flag, 16) != 1:
                direction = '+'
            else:
                direction = '-'
            mapped_coord = int(row[3])
            if AS > 0.6 * temp_fastq_length:
                result.append([mapped_coord, direction])
        return result

    def measure_length_main(self, coord1, coord2):
        """Return the distance between two coordinates of rDNA in the given read.

        Coordinate1 should be smaller than coordinate 2.

        Args:
            header (str): header
            read (str): nanopore read
            coord1 (int): coordinate in rDNA
            coord2 (int): coordinate in rDNA

        Returns:
            distances: list of distances between the given coordinates in the read
        """
        # extend items with [1] or [2] to track their coordinate.
        coord1_in_read = [i + [1] for i in _find_coordinate(read, coord1)]
        coord2_in_read = [i + [2] for i in _find_coordinate(read, coord2)]
        combined_coord = coord1_in_read + coord2_in_read

        distances = []
        flag = 0
        previous = ''
        for item in sorted(combined_coord):
            if flag == 0:
                previous = item
                flag = 1
            else:
                if previous[2] == 1:
                    if previous[1] == '+' and item[1] == '+' and item[2] == 2:
                        distances.append(item[0] - previous[0])
                        flag = 0
                    else:
                        previous = item
                        flag = 1
                else:
                    if previous[1] == '-' and item[1] == '-' and item[2] == 1:
                        distances.append(item[0] - previous[0])
                        flag = 0
                    else:
                        previous = item
                        flag = 1
        return distances


def find_rDNA_reads(fastq_path, output_name, ref):
    """Extracts rDNA reads from mixed read fastqs.

    Args:
        fastq_path (str): path of fastq
        output_name (str): the name of output fastq file, containing only rDNA reads
        ref (class): reference class
    Returns:
        None (writes rDNA fastq to output_name)
    """
    ret_str = ''
    with open(fastq_path) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            r = Fastqread(items, ref.path)
            if r.health != 0:
                ret_str += ''.join(items)
    with open(output_name, 'w') as fw:
        fw.write(ret_str)


def get_distances(fastq, coord1, coord2, save_filename='figs/length_distribution.png'):
    """Plots distance distribution.

    Args:
        fastq_path (str): path of fastq
        coord1 (int): smaller coordinate
        coord2 (int): larger coordinate
        save_filename (str): saved filename
    Returns:
        None (plots distribution)
    """
    ref_dist = coord2 - coord1
    distances = []
    with open(fastq) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            read = Fastqread(items)
            distances.append(read.measure_length_main(coord1, coord2))
    plt.hist(distances, range=(ref_dist-3000, ref_dist+3000), bins=50)
    plt.vlines([ref_dist], 0, 2000, 'black', linestyle='dashed', linewidth=0.6)
    plt.title('Length distribution {}-{}'.format(coord1, coord2))
    plt.savefig(save_filename, dpi=300)
    plt.close()


def plot_distance_differences(fastq, coord1, coord2, save_filename='figs/length_diff_distribution.png'):
    """Plot the difference between the distances of two rDNA copies.

    This function plots two histograms. One for neighboring rDNA copies and
    another for the distance difference of two randomly chosen rDNA copies.

    Args:
        fastq_path (str): path of fastq
        coord1 (int): smaller coordinate
        coord2 (int): larger coordinate
        save_filename (str): saved filename
    Returns:
        None (plots distribution)
    """
    ref_dist = coord2 - coord1
    distances = []
    with open(fastq) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            read = Fastqread(items)
            distances.append(read.measure_length_main(coord1, coord2))
    len_dis = len(distances)
    simulated_diffs = []
    for n in range(10000):
        d1 = distances[int(random.random() * len_dis)]
        d2 = distances[int(random.random() * len_dis)]
        simulated_diffs.append(abs(d1 - d2))
    actual_diffs = []
    for item in distances:
        if len(item[1]) > 1 and all(i < 30000 for i in item[1]):
            for n in range(len(item[1]) - 1):
                actual_diffs.append(abs(item[1][n] - item[1][n+1]))

    plt.subplot(2, 1, 1)
    plt.hist(actual_diffs, bins=100, range=(0, 10000))
    plt.title('Actual distribution')
    plt.subplot(2, 1, 2)
    plt.hist(simulated_diffs[:len(actual_diffs)], bins=100, range=(0, 10000))
    plt.title('Simulated distribution (' + str(len(actual_diffs)) + ')')
    plt.tight_layout()
    plt.savefig(save_filename, dpi=300)
    plt.close()


def single_fast5_to_fastqread(filename, ref):
    """Generates Fastqread object from a fast5 file. 

    Args:
        fastq_path (str): path of fast5
        ref (object): reference object
    Returns:
        read (Fastqread object)
    """
    with h5py.File(fast5file) as f:
        mod_loc = self._find_basecall_with_mod(f)
        fastq = f['/Analyses/Basecall_1D_00' + mod_loc + '/BaseCalled_template/Fastq']\
                [()].decode().strip().split('\n')
    return Fastqread(fastq, ref)


if __name__ == '__main__':
    target_id = 'fda2423a-5429-4625-a946-12f2f15bdee5'
    normal_ref = Reference('rDNA_index/humRibosomal.fa')
    with open('clive_rDNA_reads2.fastq') as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            if target_id in items[0]:
                r = Fastqread(items, normal_ref)
                break
    inv_pos = 87291
    with open('temp_files/before.fq', 'w') as fw:
        fw.write('>before\n{}\n+\n{}'.format(r.seq[inv_pos-100:inv_pos], r.quality[inv_pos-100:inv_pos]))
    with open('temp_files/after.fq', 'w') as fw:
        alen = 150
        fw.write('>after\n{}\n+\n{}'.format(r.seq[inv_pos:inv_pos+alen], r.quality[inv_pos:inv_pos+alen]))
    utilities.bwa_mapping(normal_ref.path, 'temp_files/before.fq', 'temp_files/before.sam')
    utilities.bwa_mapping(normal_ref.path, 'temp_files/after.fq', 'temp_files/after.sam')
    quit()
    print([ord(i)-33 for i in r.quality[inv_pos-20:inv_pos]])
    print([ord(i)-33 for i in r.quality[inv_pos:inv_pos+20]])
    quit()
    shutil.rmtree('temp_figs')
    os.mkdir('temp_figs')
    with open('inverted_reads.fastq') as f:
        for n, item in enumerate(itertools.zip_longest(*[iter(f)]*4)):
            r = Fastqread(item, normal_ref)
            r.plot_structure(ref, 'temp_figs')

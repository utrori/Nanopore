import random
import os
import subprocess
from matplotlib import pyplot as plt
from search_rDNA_reads import make_temp_fastq
from rDNA_structure_of_long_reads import easy_flag
import pandas as pd
import itertools


with open('rDNA_index/humRibosomal.fa') as f:
    rDNA = ''.join((line.strip() for line in f.readlines()[1:]))
FNULL = open(os.devnull, 'w')


def is_rDNAs_healthy(samdata):
    """Check if the read is completely inside rDNA or at the end.

    Args:
        samdata (list(str, str, ...)): samfile for the split mapped read
    Returns:
        healthy: 1 if so 0 if not
    """
    mapping_state = []
    for read in samdata:
        row = read.split()
        flag = int(row[1])
        healthy = 0
        if easy_flag(flag, 4) != 1:
            mapping_state.append(1)
        else:
            mapping_state.append(0)
    threshold = 0.9
    if sum(mapping_state)/len(mapping_state) > threshold:
        healthy = 1
    else:
        for i in range(1, len(mapping_state) - 100):
            if sum(mapping_state[i:])/len(mapping_state[i:]) > threshold or \
               sum(mapping_state[:-i])/len(mapping_state[:-i]) > threshold:
                healthy = 1
                return healthy
    return healthy


def make_temp_bwa_index(header, read):
    """Make bwa index from the reads.

    Args:
        header (str): header
        read (str): nanopore read
    Returns:
        None. This script just makes fasta and index by bwa.
    """
    with open('temp_index/temp_index.fasta', 'w') as fw:
        fw.write('>' + header.split()[0] + '\n' + read)
    subprocess.run('bwa index temp_index/temp_index.fasta', shell=True,
                   stdout=FNULL, stderr=subprocess.STDOUT)


def find_coordinate(read, coordinate):
    """Find exact rDNA coordinate in the given nanopore reads.

    Args:
        read (str): read sequence
        coordinate (int): target coordinate of rDNA

    Returns:
        result: list of [mapped coordinate, direction of mapping]
    """
    result = []
    temp_fastq_length = 500
    with open('coordinate_rDNA.fastq', 'w') as fw:
        fw.write('>temp\n' + rDNA[coordinate-1:coordinate+temp_fastq_length-1]
                 + '\n+\n' + 'J' * temp_fastq_length + '\n')
    # with -a option, multiple hits are more clearly shown
    subprocess.run('bwa mem -Ma -x ont2d -t '
                   '/home/yutaro/nanopore/clive/temp_index/temp_index.fasta '
                   'coordinate_rDNA.fastq > temp_sam4coord.sam', shell=True,
                   stdout=FNULL, stderr=subprocess.STDOUT)
    with open('temp_sam4coord.sam') as samf:
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


def measure_length_main(header, read, coordinate1, coordinate2):
    """Return the distance between two coordinates of rDNA in the given read.

    Coordinate1 shoudl be smaller than coordinate 2.

    Args:
        header (str): header
        read (str): nanopore read
        cooridnate1 (int): coordinate in rDNA
        cooridnate2 (int): coordinate in rDNA

    Returns:
        distances: list of distances between the given coordinates in the read
    """
    make_temp_bwa_index(header, read)
    # extend items with [1] or [2] to track their coordinate.
    coord1_in_read = [i + [1] for i in find_coordinate(read, coordinate1)]
    coord2_in_read = [i + [2] for i in find_coordinate(read, coordinate2)]
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


def measure_length_in_fastq(fastq, split_length, coord1, coord2):
    """Measure the distances between 2 coordinates in rDNA for a set of reads.

    This module takes a .fastq file that contains many multiple rDNA containing
    reads and measures the length for each of them and returns the result.
    Coordinate1 should be smaller than coordinate 2.

    Args:
        fastq (str): fastq filename
        split_length (int): split length when the read is mapped to rDNA
        coordinate1 (int): coordinate in rDNA
        coordinate2 (int): coordinate in rDNA
    Returns:
        result: list of (header, lisf of measured distance between 2 coords)
    """
    with open(fastq) as f:
        result = []
        for data in itertools.zip_longest(*[iter(f)]*4):
            header = data[0].strip()
            read = data[1].strip()
            quality = data[3].strip()
            make_temp_fastq(split_length, header, read, quality)
            subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/'
                           'clive/rDNA_index/humRibosomal.fa temp_fastq.fastq'
                           ' > temp_sam.sam', shell=True, stdout=FNULL,
                           stderr=subprocess.STDOUT)
            with open('temp_sam.sam') as samf:
                samdata = samf.readlines()[2:]
            if is_rDNAs_healthy(samdata):
                dist = measure_length_main(header, read, coord1, coord2)
                result.append((header.split()[0], dist))
    return result


def plot_distance_differences(result):
    """Plot the difference between the distances of two rDNA copies.

    This function plots two histograms. One for neighboring rDNA copies and
    another for the distance difference of two randomly chosen rDNA copies.

    Args:
        result (list(str, list)): result from measure_length_main
    Returns:
        None
    """
    distances = []
    for item in result:
        for dist in item[1]:
            if dist < 30000:
                distances.append(dist)
    len_dis = len(distances)
    simulated_diffs = []
    for n in range(10000):
        d1 = distances[int(random.random() * len_dis)]
        d2 = distances[int(random.random() * len_dis)]
        simulated_diffs.append(abs(d1 - d2))
    actual_diffs = []
    for item in result:
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
    plt.savefig('figs/noncoding_length_differences.png', dpi=300)


def plot_distances(result):
    """Plot the distances of two coordinates in rDNA for a set of reads.

    Args:
        result (list(str, list)): result from measure_length_main
    Returns:
        None
    """
    distances = []
    for item in result:
        for dist in item[1]:
            if dist < 30000:
                distances.append(dist)
    plt.hist(distances, range=(13000, 22000), bins=50)
    plt.vlines([14000], 0, 2000, 'black', linestyle='dashed', linewidth=0.6)
    plt.title('Noncoding length distribution (' + str(len(distances)) + ')')
    plt.savefig('figs/length_distribution_unstable_region.png', dpi=300)


if __name__ == '__main__':
    split_length = 200
    #result = measure_length_in_fastq('rDNA_reads.fastq', split_length, 10000, 18000)
    #pd.to_pickle(result, 'length_distribution_after_coding.pkl')
    result = pd.read_pickle('length_distribution.pkl')
    plot_distance_differences(result)
    plt.close()
    plot_distances(result)
    plt.close()
    quit()
    for item in result:
        if any(i < 14300 for i in item[1]):
            print(item)

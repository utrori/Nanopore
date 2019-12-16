import random
import os
import subprocess
from matplotlib import pyplot as plt
from search_rDNA_reads import split_sequence, make_temp_fastq
from rDNA_structure_of_long_reads import easy_flag
import pandas as pd
import itertools


with open('rDNA_index/humRibosomal.fa') as f:
    rDNA = ''.join((line.strip() for line in f.readlines()[1:]))
FNULL = open(os.devnull, 'w')


# either the read is completely inside rDNA or the read is at the end of rDNA
def is_rDNAs_healthy(samdata):
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
            if sum(mapping_state[i:])/len(mapping_state[i:]) > threshold or sum(mapping_state[:-i])/len(mapping_state[:-i]) > threshold:
                healthy = 1
                return healthy
    return healthy


def find_healthy_coord(coords, target_coord, split_length):
    margin = split_length * 5
    within_range = []
    for coord in coords:
        if (coord[1] - margin) < target_coord < (coord[1] + margin):
            within_range.append(coord)


def measure_unstable_length(samdata, split_length, read):
    coords = []
    for n, split_read in enumerate(samdata):
        row = split_read.split()
        if easy_flag(int(row[1]), 4) != 1:
            coords.append(n, int(row[3]))
        else:
            coords.append(n, 6000)
            # becuase coordinate 6000 can be ignored in our analysis, which focuses on NC region.
    find_healthy_coord(coords, 17000, split_length)


def make_temp_bwa_index(header, read):
    with open('temp_index/temp_index.fasta', 'w') as fw:
        fw.write('>' + header.split()[0] + '\n' + read)
    subprocess.run('bwa index temp_index/temp_index.fasta', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)


def find_coordinate(read, coordinate):
    result = []
    temp_fastq_length = 500
    with open('coordinate_rDNA.fastq', 'w') as fw:
        fw.write('>temp\n' + rDNA[coordinate-1:coordinate+temp_fastq_length-1] + '\n+\n' + 'J' * temp_fastq_length + '\n')
    # with -a option, multiple hits are more clearly shown
    subprocess.run('bwa mem -Ma -x ont2d -t 5 /home/yutaro/nanopore/clive/temp_index/temp_index.fasta coordinate_rDNA.fastq > temp_sam4coord.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
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


#fw.write(header.split()[0] + '\n' + read + '\n+\n' + 'J' * len(read) + '\n')
def measure_length_main(header, read, coordinate1, coordinate2):
    # coord1 should be smaller than coord2
    make_temp_bwa_index(header, read)
    coord1_in_read = [item + [1] for item in find_coordinate(read, coordinate1)]
    coord2_in_read = [item + [2] for item in find_coordinate(read, coordinate2)]
    combined_coord = coord1_in_read + coord2_in_read
    sorted_combined = sorted(combined_coord)
    # automatically sorted by the first element!

    distances = []
    flag = 0
    previous = ''
    for item in sorted_combined:
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


def measure_length_in_fastq(fastq, split_length, coordinate1, cooridnate2):
    with open(fastq) as f:
        result = []
        for data in itertools.zip_longest(*[iter(f)]*4):
            header = data[0].strip()
            read = data[1].strip()
            quality = data[3].strip()
            make_temp_fastq(split_length, header, read, quality)
            subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_fastq.fastq > temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            with open('temp_sam.sam') as samf:
                samdata = samf.readlines()[2:]
            if is_rDNAs_healthy(samdata):
                result.append((header.split()[0], measure_length_main(header, read, coordinate1, cooridnate2)))
    return result


def plot_distance_differences(result):
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
    plt.title('Simulated distribution')
    plt.tight_layout()
    
    plt.savefig('distance_differences.png')


def plot_distances(result):
    distances = []
    for item in result:
        for dist in item[1]:
            if dist < 30000:
                distances.append(dist)
    plt.hist(distances, range=(12000, 25000), bins=50)
    plt.vlines([14000], 0, 2000, 'black', linestyle='dashed', linewidth=0.6)
    plt.savefig('length_distribution_unstable_region.png', dpi=300)


if __name__ == '__main__':
    split_length = 200
    #result = measure_length_in_fastq('rDNA_reads.fastq', split_length, 18000, 32000)
    #pd.to_pickle(result, 'length_distribution.pkl')
    result = pd.read_pickle('length_distribution.pkl')
    plot_distance_differences(result)
    quit()
    for item in result:
        if any(i > 22000 for i in item[1]):
            print(item)

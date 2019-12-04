# -*- coding: utf-8 -*-
# This script focuses on longer reads and determine the direction and the distance between repeats.

from search_rDNA_reads import split_sequence, make_temp_fastq
import subprocess
import os
import pandas as pd


def easy_flag(flag, base):
    return int(flag % (2 * base) / base)


def analyze_split_reads(samfile, split_length):
    rDNA_len = 42999
    rDNA_TR_len = 13314
    flags = []
    with open(samfile) as sf:
        samdata = sf.readlines()[2:]
        # in_TR is a flag used to keep track of the where the split read is (is it inside TR o not). Not to confuse with flag (SAM header flag).
        in_TR = 0
        # TRs is a colletcion of read sets. Each read set is composed of reads that are mapped to each TR.
        TRs = []
        TR = []
        prev_direction = ''
        for n, item in enumerate(samdata):
            row = item.split()
            flag = int(row[1])
            if easy_flag(flag, 4) == 1:
            # reads that do not map to rDNA
                continue
            if easy_flag(flag, 256) == 1:
            # reads that are mapped to multiple loci
            # temoporarlily disregard them because most of them should be in the non-transcribed region (repeats in the NTR).
                continue
            if easy_flag(flag, 16) == 1:
                direction = '-'
            else:
                direction = '+'
            if direction == '+':
                position = [int(row[3]), (int(row[3]) + split_length) % rDNA_len]
            else:
                position = [(int(row[3]) - split_length) % rDNA_len, int(row[3])]
            if 0 < position[0] < rDNA_TR_len and 0 < position[1] < rDNA_TR_len:
                # when the reads change its direction
                if prev_direction != '' and prev_direction != direction and len(TR) > 10:
                    TRs.append(TR)
                    TR = []
                in_TR = 1
                # TR definition
                TR.append((n, position[0], position[1], direction))
            else:
                in_TR = 0
                if len(TR) > 10:
                    TRs.append(TR)
                TR = []
            prev_direction = direction
    # end position
    TRs.append([(n, 0, 0, '*')])
    return TRs


def visualize_TRs(TRs, split_length):
    size_of_bin = 1000
    prev_end = 0
    TR_lens = []
    dists = []
    directions = []
    prev_direction = ''
    output = ''
    if TRs == []:
        return output
    for TR in TRs[:-1]:
        # start and end are the numbers of split
        start = TR[0][0]
        direction = TR[0][3]
        directions.append(direction)
        end = TR[-1][0]
        TR_lens.append((end - start) * split_length)
        dists.append((start - prev_end) * split_length)
        prev_end = end
        prev_direction = direction
    for TR_len, dist, direction in zip(TR_lens, dists, directions):
        if direction == '+':
            output += '-' * int(dist / size_of_bin) + '￫' * int(TR_len/ size_of_bin)
        else:
            output += '-' * int(dist / size_of_bin) + '￩' * int(TR_len/ size_of_bin)
    dist = (TRs[-1][0][0] - prev_end) * split_length
    output += '-' * int(dist / size_of_bin)
    return output


def analyze_direction_distance(TRs, split_length):
    result = []
    for n in range(len(TRs) - 2):
        if TRs[n][0][3] == TRs[n+1][0][3]:
            direction = '='
        else:
            direction = '!'
        distance = split_length * (TRs[n+1][-1][0] - TRs[n][-1][0])
        result.append((direction, distance))
    return result


def count_direction_dist(set_of_TRs, split_length):
    summary = []
    for TRs in set_of_TRs:
        summary.extend(analyze_direction_distance(TRs[1], split_length))
    same = 0
    same_distances = []
    inverse = 0
    inverse_distances = []
    for item in summary:
        dist = item[1]
        if item[0] == '=':
            same += 1
            same_distances.append(dist)
        else:
            inverse += 1
            inverse_distances.append(dist)
    print((same, inverse))
    print(inverse_distances)


def analyze_all_fastq_file(fastq, split_length, length_cutoff):
    count = 0
    set_of_TRs = []
    with open(fastq) as f:
        for n, line in enumerate(f):
            print(n)
            if n % 4 == 0:
                header = line.strip()
            if n % 4 == 1:
                read = line.strip()
            if n % 4 == 3:
                quality = line.strip()
                if len(quality) < length_cutoff:
                    continue
                else:
                    count += 1
                    make_temp_fastq(split_length, header, read, quality)
                    FNULL = open(os.devnull, 'w')
                    subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_fastq.fastq > temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    TRs = analyze_split_reads('temp_sam.sam', split_length)
                    set_of_TRs.append([header, TRs])
    return set_of_TRs


def print_TRs(TRs):
    for TR in TRs:
        print(TR)


def find_TRs_by_header(set_of_TRs, header_id):
    for item in set_of_TRs:
        if header_id in item[0]:
            return item[1]


def analyze_fastq_by_header(fastq, header_id, split_length):
    with open(fastq) as f:
        fourd = 0
        for line in f:
            if header_id in line:
                header = line
                read = f.readline()
                f.readline()
                quality = f.readline()
                make_temp_fastq(split_length, header, read, quality)
                FNULL = open(os.devnull, 'w')
                subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_fastq.fastq > single_split_mapped.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                break


if __name__ == '__main__':
    split_length = 200
    length_cutoff = 80000
    # this set_of_TRs is a set of [header, TRs]
    # TRs is a list of split reads that are mapped to the Trancribed Region. The number of the last split reads is added to the end of it. (see line 59)
    #set_of_TRs = analyze_all_fastq_file('rDNA_reads.fastq', split_length, length_cutoff)
    #pd.to_pickle(set_of_TRs, 'set_of_TRs.pkl')
    set_of_TRs = pd.read_pickle('set_of_TRs.pkl')
    print_TRs(set_of_TRs[0][1])
    count_direction_dist(set_of_TRs, split_length)
    analyze_fastq_by_header('rDNA_reads.fastq', '@fda2423a-5429-4625-a946-12f2f15bdee5', split_length)
    print_TRs(find_TRs_by_header(set_of_TRs, '@fda2423a-5429-4625-a946-12f2f15bdee5'))
    with open('reads_visualized.txt', 'w') as fw:
        for n, TRs in enumerate(set_of_TRs):
            fw.write(str(n) + '\t' + TRs[0].split()[0] + '\n' +  visualize_TRs(TRs[1], split_length) + '\n')

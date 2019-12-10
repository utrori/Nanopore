from search_rDNA_reads import split_sequence, make_temp_fastq
from rDNA_structure_of_long_reads import easy_flag
import os
import subprocess
import itertools
import pandas as pd


def find_end_reads(samfile, split_length):
    with open(samfile) as f:
        samdata = f.readlines()[2:]
    mapping_state = []
    threshold = 0.2
    for item in samdata:
        row = item.split()
        flag = int(row[1])
        if easy_flag(flag, 4) != 1:
            if easy_flag(flag, 16) != 1:
                mapping_state.append(1)
            else:
                mapping_state.append(-1)
            # 1 for + direction and -1 for - direction
        else:
            mapping_state.append(0)

    offset = 20
    # you have to ignore the end of reads because it can cause very high non-mapping rate by chance
    candidates = []
    for n in range(offset, len(mapping_state) - offset):
        left = [abs(i) for i in mapping_state[:n]]
        right = [abs(i) for i in mapping_state[n:]]
        left_rate = sum(left) / len(left)
        right_rate = sum(right) / len(right)
        if left_rate > (1 - threshold) and right_rate < threshold:
            candidates.append((n, left_rate, right_rate))
        elif left_rate < threshold and right_rate > (1 - threshold):
            candidates.append((n, left_rate, right_rate))
    if candidates == []:
        return 0
    else:
        provisional_best_diff = 0
        provisional_best_pos = 0
        for n, lr, rr in candidates:
            diff = abs(lr - rr)
            if diff > provisional_best_diff:
                provisional_best_diff = diff
                provisional_best_pos = n
                provisional_lr = lr
                provisional_rr = rr

        lr = provisional_lr
        rr = provisional_rr
        n = provisional_best_pos
        if lr > rr:
            side_of_non_rDNA = 'right'
            near_boundary_directions = [i for i in mapping_state[n-20:n]]
        else:
            side_of_non_rDNA = 'left'
            near_boundary_directions = [i for i in mapping_state[n:n+20]]
        # average direction of 20 reads at the boudary region
        if sum(near_boundary_directions) > 0:
            direction = '+'
        else:
            direction = '-'
        return (split_length * n, direction, side_of_non_rDNA)


def make_fastq_for_boundary(fastq, boundary):
    half_length = 1000
    output = fastq[0] + '\n'
    output += fastq[1][boundary - half_length:boundary + half_length] + '\n'
    output += '+' + '\n'
    output += fastq[2][boundary - half_length:boundary + half_length]
    with open('temp_boundary.fastq', 'w') as fw:
        fw.write(output)


if __name__ == '__main__':
    fastq = 'rDNA_reads.fastq'
    FNULL = open(os.devnull, 'w')
    split_length = 200

    with open('rDNA_reads.fastq') as f:
        boundaries = []
        for n, each_fastq in enumerate(itertools.zip_longest(*[iter(f)]*4)):
            print(n)
            header = each_fastq[0]
            read = each_fastq[1]
            quality = each_fastq[3]
            make_temp_fastq(split_length, header, read, quality)
            subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_fastq.fastq > temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

            boundary = find_end_reads('temp_sam.sam', split_length)
            if boundary != 0:
                make_fastq_for_boundary((header, read, quality), boundary[0])
                subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_boundary.fastq > temp_boundary.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                with open('temp_boundary.sam') as samf:
                    samdata = samf.readlines()[2:]
                cigar = samdata[0].split()[5]
                # only when the CIGAR sequence is simple
                dist_to_boundary = ''
                if boundary[2] == 'right':
                    temp_dist = cigar.split('M')[-1]
                    if 'S' in temp_dist:
                        dist_to_boundary = int(temp_dist[:-1])
                        true_boundary = boundary[0] + 1000 - dist_to_boundary
                # non rDNA reads are expressed as read[true_boundary:]
                else:
                    temp_dist = cigar.split('S')[0]
                    if temp_dist.isdigit():
                        dist_to_boundary = int(temp_dist)
                        true_boundary = boundary[0] -1000 + dist_to_boundary
                # non rDNA reads are expressed as read[:true_boundary]
                if dist_to_boundary == '':
                    continue
                else:
                    boundaries.append((header.split()[0], true_boundary, boundary[1], boundary[2]))
                    # header, boundary_position, direction_of_rDNA, side_of_non_rDNA
            pd.to_pickle(boundaries, 'rDNA_boundaries.pkl')

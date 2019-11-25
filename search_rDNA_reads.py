import subprocess
import matplotlib.pyplot as plt
import os


def check_rDNA_reads_in_sam(filename):
    with open(filename) as f:
        f.readline() 
        f.readline() 
        # this alignment score is normalized by read length
        alignment_scores = []
        for line in f:
            row = line.split()
            # the second row of SAM file header is mapping status. 256 means the mapping is no primary and shorter split mapping.
            # the third row is contig the read is mapped on. If the read is not mapped to any part of the index sequence, it shows *.
            if row[2] !='*' and int(int(int(row[1]) % 512) / 256) != 1:
                read_length = len(row[9])
                # when the read_length is too short, it can produce very high score
                if read_length > 500:
                    alignment_scores.append(float(line.split('AS:i:')[1].split()[0]) / read_length)
    return alignment_scores


def split_sequence(sequence, split_length):
    split_seq = []
    for i in range(int(len(sequence) / split_length + 1)):
        split_seq.append(sequence[i*split_length:(i+1)*split_length])
    return split_seq


count = 0
# This process split each read into the size split_length and map them to rDNA consensus sequence.
# The header for each split sequence is "the original header"_number.
split_length = 4000
rDNA_containing_reads = []
alignment_scores = []
with open('/home/yutaro/data/cliveome/fastq_runid_0ee57c1a265c2f494821757929f1af60fe060a43_0.fastq') as f:
    for n, line in enumerate(f):
        if n % 4 == 0:
            print(int(n/4))
            header = line.strip()
        if n % 4 == 1:
            read = line.strip()
            split_reads = split_sequence(read, split_length)
        if n % 4 == 3:
            quality = line.strip()
            split_qualities = split_sequence(quality, split_length)
            # here write temporary .fastq made by splitting each read. 
            with open('temp_fastq.fastq', 'w') as fw:
                for i in range(len(split_reads)):
                    split_header = [header.split()[0], ' '.join(header.split()[1:])]
                    fw.write(split_header[0] + '_' + str(i+1) + ' ' + split_header[1] + '\n' + split_reads[i] + '\n+\n' + split_qualities[i] + '\n')
            # perform bwa for the created fastq
            # you can avoid printing the output of subprocess by directing the stdout to devnull
            FNULL = open(os.devnull, 'w')
            subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_fastq.fastq > temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            alignment_scores = check_rDNA_reads_in_sam('temp_sam.sam')
            if any(AS > 0.5 for AS in alignment_scores):
                print(alignment_scores)
                with open('rDNA_containing_reads.fastq', 'a') as fw:
                    fw.write(header + '\n' + read + '\n+\n' + quality + '\n')

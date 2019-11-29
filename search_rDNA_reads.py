import subprocess
import sys
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
            # the second row of SAM file header is mapping status. 256 means the mapping is shorter split mapping, not primary.
            # the third row is contig on which the read is mapped. If the read is not mapped to any part of the index sequence, it shows *.
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


# creates a temporary .fastq made by splitting each read.
# The header for each split sequence is "the original header"_number.
def make_temp_fastq(split_length, header, read, quality, tempfilename='temp_fastq.fastq'):
    split_reads = split_sequence(read, split_length)
    split_qualities = split_sequence(quality, split_length)
    with open('temp_fastq.fastq', 'w') as fw:
        for i in range(len(split_reads)):
            split_header = [header.split()[0], ' '.join(header.split()[1:])]
            fw.write(split_header[0] + '_' + str(i+1) + ' ' + split_header[1] + '\n' + split_reads[i] + '\n+\n' + split_qualities[i] + '\n')


if __name__ == '__main__':
    # This process split each read into the size of split_length and map them to rDNA consensus sequence.
    split_length = 4000
    #fastq for test
    #input_fastq = '/home/yutaro/data/cliveome/fastq_runid_0ee57c1a265c2f494821757929f1af60fe060a43_0.fastq'
    input_fastq = sys.argv[1]
    with open(input_fastq) as f:
        for n, line in enumerate(f):
            if n % 4 == 0:
                header = line.strip()
            if n % 4 == 1:
                read = line.strip()
            if n % 4 == 3:
                quality = line.strip()
                make_temp_fastq(header, read, quality)

                # perform bwa for the created fastq
                # you can avoid printing the output of subprocess by directing the stdout to devnull
                FNULL = open(os.devnull, 'w')
                subprocess.run('bwa mem -M -x ont2d -t 5 /home/yutaro/nanopore/clive/rDNA_index/humRibosomal.fa temp_fastq.fastq > temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                alignment_scores = check_rDNA_reads_in_sam('temp_sam.sam')
                if any(AS > 0.5 for AS in alignment_scores):
                    with open('rDNA_containing_reads.fastq', 'a') as fw:
                        fw.write(header + '\n' + read + '\n+\n' + quality + '\n')

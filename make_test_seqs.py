# This script makes an artificial .fastq file to submit to bwa mem.
# Here I use mouse ribosomal rRNA consensus sequence.

import random


def simulating_read_errors(seq):
    mutated_seq = ''
    for n in seq:
        if random.random() > 0.80:
            mutated_seq += 'A'
        else:
            mutated_seq += n
    return mutated_seq


read_length = 4000
filename_mus = '/home/yutaro/data/Mus_musculus_NCBI_GRCm38/Mus_musculus/NCBI/GRCm38/Sequence/AbundantSequences/musRibosomal.fa'
filename_homo = '/home/yutaro/data/Homo_sapiens_NCBI_GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/AbundantSequences/humRibosomal.fa'
with open(filename_homo) as f:
    f.readline()
    rDNA_seq = ''
    for line in f:
        rDNA_seq += line.strip()
    sample_seqs = []
    for i in range(10):
        sample_seqs.append(rDNA_seq[i*read_length:(i+1)*read_length])

headers = []
qualities = []

# The quality track is obtained from a cliveome file.
with open('/home/yutaro/data/cliveome/fastq_runid_0ee57c1a265c2f494821757929f1af60fe060a43_0.fastq') as f:
    for n, line in enumerate(f):
        if n % 4 == 0:
            headers.append(line.strip())
        if n % 4 == 3:
            if len(line) > read_length:
                qualities.append(line[0:read_length])
        if n > 200:
            break

with open('test_data.fastq', 'w') as fw:
    for i in range(10):
        fw.write(headers[i] + '\n' + sample_seqs[i] + '\n+\n' + qualities[i] + '\n')

with open('test_data_mutated.fastq', 'w') as fw:
    for i in range(10):
        fw.write(headers[i] + '\n' + simulating_read_errors(sample_seqs[i]) + '\n+\n' + qualities[i] + '\n')



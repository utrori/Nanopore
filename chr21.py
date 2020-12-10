import subprocess
import os
import itertools

chr21 = ""

with open('/home/yutaro/data/Homo_sapiens_NCBI_GRCh38/Homo_sapiens/NCBI/GRCh38/Sequence/Chromosomes/chr21.fa') as f:
    f.readline()
    for line in f:
        chr21 += line.strip()

print(len(chr21))
print(chr21[:50])

with open('split_chr21.fastq', 'w') as fw:
    for n in range(int(len(chr21) / 100)):
        fw.write('@chr21_' + str(n) + '\n')
        fw.write(chr21[n * 100:(n+1) * 100] + '\n+\n' + 'J' * 100 + '\n')

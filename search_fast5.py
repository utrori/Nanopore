import itertools


read_ids = []
with open('rDNA_containing_reads.fastq') as f:
    for n, item in enumerate(itertools.zip_longest(*[iter(f)]*4)):
        read_ids.append(item[0].split()[0][1:])

read_ids2 = []
with open('fastq2fast5.txt') as f:
    for line in f:
        read_ids2.append(line.strip().split()[1])

for i in read_ids:
    if i not in read_ids2:
        print(i)
quit()
retstr = ''
with open('sequencing_summary.txt') as f:
    for n, line in enumerate(f):
        print(len(read_ids))
        row = line.split()
        if row[1] in read_ids:
            retstr += row[0] + '\t' + row[1] + '\n'
            read_ids.remove(row[1])

with open('fastq2fast5.txt', 'w') as fw:
    fw.write(retstr)

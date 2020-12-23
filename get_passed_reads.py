import itertools

ids_in_summary = []
with open('sequencing_summary.txt') as f:
    f.readline()
    for line in f:
        row = line.split()
        ids_in_summary.append(row[1])

ret_str = ''
with open('clive_rDNA_reads2.fastq') as f:
    for items in itertools.zip_longest(*[iter(f)]*4):
        if items[0].split()[0][1:] in ids_in_summary:
            ret_str += ''.join(items)
with open('clive_rDNA_reads3.fastq', 'w') as fw:
    fw.write(ret_str)

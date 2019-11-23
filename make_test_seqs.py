# This script makes an artificial .fastq file to submit to bwa mem.
# Here I use mouse ribosomal rRNA consensus sequence.

read_length = 1000
with open('/home/yutaro/data/Mus_musculus_NCBI_GRCm38/Mus_musculus/NCBI/GRCm38/Sequence/AbundantSequences/musRibosomal.fa') as f:
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

with open('test_data_mus.fastq', 'w') as fw:
    for i in range(10):
        fw.write(headers[i] + '\n' + sample_seqs[i] + '\n+\n' + qualities[i] + '\n')

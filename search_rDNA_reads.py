import subprocess


def check_rDNA_reads_in_sam(filename):
    with open(filename) as f:
       f.readline() 
       f.readline() 


# This process split each read into the size split_length and map them to rDNA consensus sequence.
# The header for each split sequence is "the original header"_number.
split_length = 4000
with open('/home/yutaro/data/cliveome/fastq_runid_0ee57c1a265c2f494821757929f1af60fe060a43_0.fastq') as f:
    for n, line in enumerate(f):
        if n % 4 == 0:
            header = line.strip()
        if n % 4 == 1:
            read = line.strip()
            split_reads = []
            for i in range(int(len(read) / split_length)):
                split_reads.append(read[i*split_length:(i+1)*split_length])
        if n % 4 == 3:
            quality = line.strip()
            split_qualities = []
            for i in range(int(len(read) / split_length)):
                split_qualities.append(quality[i*split_length:(i+1)*split_length])
            # here write temporary .fastq made by splitting each read. 
            with open('temp_fastq.fastq', 'w') as fw:
                for i in range(len(split_reads)):
                    split_header = [header.split()[0], ' '.join(header.split()[1:])]
                    fw.write(split_header[0] + '_' + str(i+1) + ' ' + split_header[1] + '\n' + split_reads[i] + '\n+\n' + split_qualities[i] + '\n')
            # perform bwa for the created fastq
            subprocess.run('bwa mem -x ont2d -t 5 ./rDNA_index/humRibosomal.fa temp_fastq.fastq > temp_sam.sam')

# This script takes a .fastq file and calculate the read lengths of reads contained in the file.

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


filename = '/home/yutaro/data/cliveome/fastq_runid_0ee57c1a265c2f494821757929f1af60fe060a43_100.fastq'

with open(filename) as f:
    read_lengths = []
    for n, line in enumerate(f):
        if n % 4 == 1:
            read_lengths.append(len(line))


# total number of bases in the file.
print(sum((int(i) for i in read_lengths)))

# figure of read length distribution.
plt.hist(read_lengths, bins=20)
plt.show()

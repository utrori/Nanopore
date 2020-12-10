from Bio.Seq import Seq
import subprocess
import shutil
import numpy as np
import matplotlib.pyplot as plt
import glob
import itertools
import subprocess
import os
from main_analysis_fast5 import Fast5read

FNULL = open(os.devnull, 'w')

def find_2dlike_reads():
    fastqs = glob.glob('/mnt/data/cliveome/*.fastq')
    n = 0
    m = 0
    print('read_id\tread_length\tAS\tfrequency_of_2dlike')
    for f in fastqs:
        with open(f) as fh:
            for item in itertools.zip_longest(*[iter(fh)]*4):
                n += 1
                read_id = item[0].split()[0]
                seq = item[1].strip()
                quality = item[3].strip()
                boundary = int(len(seq)/2)
                former = seq[:boundary]
                latter = seq[boundary:]
                formerq = quality[:boundary]
                latterq = quality[boundary:]
                with open('temp_files/temp_ref.fa', 'w') as fw:
                    fw.write('>temp\n' + former)
                with open('temp_files/temp_fastq.fastq', 'w') as fw:
                    fw.write('@temp_seq\n' + latter + '\n+\n' + latterq)
                subprocess.run('~/Softwares/minimap2/minimap2 -a ' + 'temp_files/temp_ref.fa temp_files/temp_fastq.fastq > temp_files/temp_sam.sam', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                with open('temp_files/temp_sam.sam') as f:
                    result = f.readlines()[2]
                    row = result.split()
                    if row[2] != '*':
                        AS = int(result.split('AS:i:')[1].split()[0])
                        if AS/len(seq) > 0.5:
                            m += 1
                            print(read_id + '\t' + str(len(seq)) + '\t' + str(AS) + '\t' + str(m/n))
                if n > 100000:
                    quit()

ref = 'rDNA_index/humRibosomal.fa'
temp_figs = glob.glob('inverted_reads_figs/*.png')
inverse_read_ids = [i.split('/')[1].split('.')[0] for i in temp_figs]
n = 0
if os.path.exists('temp_figs'):
    shutil.rmtree('temp_figs')
os.mkdir('temp_figs')
for n, read_id in enumerate(inverse_read_ids):
    read_path = 'guppy_with_mod/workspace/' + read_id +'.fast5'
    if os.path.exists(read_path):
        if '44f76' not in read_id:
            continue
        print(read_id)
        read = Fast5read(read_path, ref)
        window = 5000
        for n, item in enumerate(read.sam_summary):
            print(item)
            if item[1] == '+':
                after_changeing_direc = n
                break
        central = after_changeing_direc * 300
        with open('temp_files/temp_fastq.fastq', 'w') as fw:
            fw.write('>temp\n' + read.seq[central -500:central+500] + '\n+\n' + read.quality[central-500:central+500])
        subprocess.run('bwa mem -M -x ont2d -t 6 ' + ref + ' temp_files/temp_fastq.fastq > temp_files/temp_sam.sam', shell=True)

        """
        x = []
        y = []
        for i in range(0, len(read.raw), window):
            x.append(i)
            y.append(np.std(read.raw[i:i+window]))
        """
        """
        for n in range(int(len(read.raw)/window)):
            plt.plot([i for i in range(window)], read.raw[n*window:(n+1)*window], linewidth=0.5)
            plt.savefig('temp_figs/' + read_id + '_' + str(n) +'.png', dpi=500)
            plt.close()
        """
        quit()

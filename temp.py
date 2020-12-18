import subprocess
import shutil
import glob
import itertools

read_ids = []
figs = glob.glob('inverted_reads_figs/*.png')
read_ids = [i.split('/')[1].split('.')[0] for i in figs]
ret_str = ''

with open('chopped_inverted.fastq') as f:
    for item in itertools.zip_longest(*[iter(f)]*4):
        print(item[0].split()[0])
        print(len(item[1]))
quit()
with open('clive_rDNA_reads.fastq') as f:
    for item in itertools.zip_longest(*[iter(f)]*4):
        if item[0].split()[0][1:] in read_ids:
            ret_str += ''.join(item)
with open('inverted_reads.fastq', 'w') as fw:
    fw.write(ret_str)
quit()
print(read_ids)
print(len(read_ids))
for rid in read_ids:
    try:
        shutil.copy('guppy_with_mod/workspace/' + rid + '.fast5', 'temp_fast5s/')
    except:
        continue
quit()
summ = 0
fs = glob.glob('/mnt/data/cliveome/*.fastq')
for f in fs:
    s = subprocess.check_output('wc -l ' + f, shell=True)
    summ += int(s.decode().split()[0])/4
print(summ)

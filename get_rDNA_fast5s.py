import pandas as pd
import itertools
import ont_fast5_api
import os
import shutil
import subprocess
import glob

ret_str = ''
read_ids = []
with open('clive_rDNA_reads2.fastq') as f:
    for item in itertools.zip_longest(*[iter(f)]*4):
        read_ids.append(item[0].split()[0][1:])

fqs = glob.glob('/mnt/data/cliveome/*.fastq')
for fq in fqs:
    print(fq)
    with open(fq) as f:
        for item in itertools.zip_longest(*[iter(f)]*4):
            read_id = item[0].split()[0][1:]
            if read_id in read_ids:
                ret_str += read_id + '\t' + fq + '\n'
with open('rDNA2fq_filename.txt', 'w') as fw:
    fw.write(ret_str)
quit()

read_id2file = pd.read_pickle('read_id2file.pkl')
print(len(read_id2file))
print(len(read_ids))
n = 0
for rid in read_ids:
    try:
        read_id2file[rid]
    except:
        n += 1
print(n)
quit()
with open('sequencing_summary.txt') as f:
    f.readline()
    for line in f:
        row = line.split()
        read_id2file[row[1]] = row[0]

pd.to_pickle(read_id2file, 'read_id2file.pkl')


id2file = pd.read_pickle('id2filename.pkl')

files = []
for read_id in read_ids:
    files.append(id2file[read_id].split('runid_')[1].split('.')[0])
files = set(files)

for file_id in files:
    print(file_id)
    for fast5 in fast5_list:
        if file_id in fast5:
            aws_com = 'aws s3 cp --no-sign-request \"s3://ont-hg1b/PromethION ' + fast5 + '\" fast5s/temp.fast5'
            subprocess.run(aws_com, shell=True)
            subprocess.run('multi_to_single_fast5 -i fast5s/temp.fast5 -s fast5s/ -t 6', shell=True)
            for filename in os.listdir('fast5s/0/'):
                if '@' + filename.split('.')[0] not in read_ids:
                    os.remove('fast5s/0/' + filename)
                else:
                    shutil.move('fast5s/0/' + filename, 'fast5s/' + filename)

            print(fast5, file_id)
            break

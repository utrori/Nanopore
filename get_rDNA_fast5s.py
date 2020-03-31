import pandas as pd
import itertools
import ont_fast5_api
import os
import shutil
import subprocess

read_ids = []
with open('reads_visualized.txt') as f:
    for item in itertools.zip_longest(*[iter(f)]*2):
        read_ids.append(item[0].split()[1].strip())

fast5_list = []
with open('fast5.listing.txt') as f:
    for line in f:
        if 'pass' in line:
            fast5_list.append(line.strip().split()[-1])

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

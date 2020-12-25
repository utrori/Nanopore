import glob
import collections
import os
import subprocess
import pandas as pd
import itertools
import matplotlib.pyplot as plt

class Clive_infos(object):
    def __init__(self):
        self.rDNA_read_ids = []
        self.fast5_pass_fns = []
        self.fast5_fail_fns = []
        self.fastq_fns = []
        self.downloaded_fastqs = []

    def get_rDNA_ids(self, fastq='clive_rDNA_reads3.fastq'):
        with open(fastq) as f:
            for items in itertools.zip_longest(*[iter(f)]*4):
                self.rDNA_read_ids.append(items[0].split()[0][1:])
        return self.rDNA_read_ids

    def get_fastq_from_listing(self):
        with open('fastq.listing.txt') as f:
            for line in f:
                if '.fastq' in line:
                    self.fastq_fns.append(line.strip().split()[-1])
        return self.fastq_fns

    def get_fastq_ids_from_listing(self):
        fastq_ids = []
        with open('fastq.listing.txt') as f:
            for line in f:
                if '.fastq' in line:
                    self.fastq_ids.append(line.strip().split('.run_id_')[-1])
        return self.fastq_ids

    def get_fast5_ids_from_listing(self):
        passed_fast5_ids = []
        failed_fast5_ids = []
        with open('fast5.listing.txt') as f:
            for line in f:
                f5id = '_'.join(line.strip().split('_')[-2:])
                if 'fast5_pass' in line:
                    passed_fast5_ids.append(f5id)
                else:
                    failed_fast5_ids.append(f5id)
        return passed_fast5_ids, failed_fast5_ids

    def get_fast5s_from_listing(self):
        with open('fast5.listing.txt') as f:
            for line in f:
                if 'fast5_pass' in line:
                    self.fast5_pass_fns.append('PromethION ' + line.strip().split()[-1])
                elif 'fast5_fail' in line:
                    self.fast5_fail_fns.append('PromethION ' + line.strip().split()[-1])
        return self.fast5_pass_fns, self.fast5_fail_fns

    def get_fastqs_from_downloaded(self):
        files = glob.glob('/mnt/data/cliveome/*.fastq')
        for f in files:
            self.downloaded_fastqs.append(f.split('/')[-1])
        return self.downloaded_fastqs

    def read_id_to_fast5(self):
        pickle_fn = 'rid2f5.pkl'
        pickle_fn2 = 'f52rid.pkl'
        if os.path.exists(pickle_fn) and os.path.exists(pickle_fn2):
            return pd.read_pickle(pickle_fn), pd.read_pickle(pickle_fn2)
        else:
            f52rid = collections.defaultdict(list)
            rid2f5 = {}
            with open('sequencing_summary.txt') as f:
                f.readline()
                for line in f:
                    row = line.split()
                    rid = row[1]
                    fast5 = row[0]
                    rid2f5[rid] = fast5
                    f52rid[fast5].append(rid)
            pd.to_pickle(rid2f5, pickle_fn)
            pd.to_pickle(f52rid, pickle_fn2)
            return rid2f5, f52rid

    def get_summary_infos(self, file_path):
        read_id2filter = {}
        read_id2length = {}
        with open(file_path) as f:
            defs = f.readline().split()
            for n in range(len(defs)):
                if 'read_id' in defs[n]:
                    read_id_pos = n
                if 'passes_filtering' in defs[n]:
                    filtering_pos = n
                if 'sequence_length_template' in defs[n]:
                    seq_len_pos = n
            for line in f:
                row = line.split()
                read_id2filter[row[read_id_pos]] = row[filtering_pos]
                read_id2length[row[read_id_pos]] = row[seq_len_pos]
        return read_id2filter, read_id2length


infos = Clive_infos()
rid2fil, rid2len = infos.get_summary_infos('sequencing_summary.txt')
pngs = glob.glob('/home/yutaro/nanopore/clive/temp_figs/*.png')
invs = [p.split('/')[-1].split('.')[0] for p in pngs]
for i in rid2fil:
    if i in invs:
        subprocess.run('grep ' + i + ' sequencing_summary.txt', shell=True)
quit()

for rid in rid2fil:
    if rid in rids3:
        print(rid)
clive_rids = infos.get_rDNA_ids

with open('sequencing_summary.txt') as f:
    for line in f:
        if 'PAD51521_1ef3eda31ddeb61dcbd8277041b49fe23a3c18df_243.fast5' in line:
            row = line.split()
            if row[7] !='TRUE':
                print(line)
quit()

quit()
ts = []
fs = []
for rid in rids:
    if r2f[rid] == 'TRUE':
        ts.append(rid)
    else:
        fs.append(rid)
print(len(ts))
print(len(fs))
quit()
f5_ids_from_listing = infos.get_fast5_ids_from_listing()
n = 0
for i in f5_ids_from_listing[0]:
    if i in f5_ids_from_listing[1]:
        n += 1
print(n)

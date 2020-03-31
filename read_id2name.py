import os
import pandas as pd
import itertools
import collections

fn2id = collections.defaultdict(list)
id2fn = {}
path = '/home/yutaro/data/cliveome/'
for fn in os.listdir(path):
    if '.fastq' in fn:
        with open(path + fn) as f:
            for read in itertools.zip_longest(*[iter(f)]*4):
                read_id = read[0].split()[0]
                fn2id[fn].append(read_id)
                id2fn[read_id] = fn

pd.to_pickle(fn2id, 'filename2id.pkl')
pd.to_pickle(id2fn, 'id2filename.pkl')

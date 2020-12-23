import h5py
import matplotlib.pyplot as plt
import utilities
import numpy as np
import matplotlib.cm as cm
import glob
import os
import shutil

class Fast5read(object):
    def _find_basecall_with_mod(self, h5file):
        for i in range(3):
            bc_template = h5file['Analyses/Basecall_1D_00' + str(i) +
                                 '/BaseCalled_template']
            if 'ModBaseProbs' in bc_template.keys():
                return str(i)

    def __init__(self, filename, ref):
        with h5py.File(filename) as f:
            read_name = list(f['Raw/Reads'].keys())[0]
            self.read_id = f['Raw/Reads/' + read_name].attrs['read_id'].decode()
            self.raw_duration = f['Raw/Reads/' + read_name].attrs['duration']
            self.raw = f['Raw/Reads/'+ read_name + '/Signal'][:]
            mod_loc = self._find_basecall_with_mod(f)
            basecalled_template = f['/Analyses/Basecall_1D_00' + mod_loc + '/BaseCalled_template']
            self.raw_start = int(f['Analyses/Segmentation_00' + mod_loc + '/Summary/segmentation'].attrs['first_sample_template'])
            self.raw_step = int(f['Analyses/Basecall_1D_00' + mod_loc + '/Summary/basecall_1d_template'].attrs['block_stride'])
            fastq = basecalled_template['Fastq'][()].decode().strip().split('\n')
            self.seq = fastq[1]
            self.quality = fastq[3]
            #fastq length! not raw signal length
            self.length = len(self.seq)
            self.guppy_mbt = basecalled_template['ModBaseProbs'][:]
            self.guppy_trace = basecalled_template['Trace'][:]
            self.guppy_move = basecalled_template['Move'][:]
            self.split_length = 300
            self.sam_summary = np.array(utilities.split_mapping_and_sam_analysis(self.split_length, 'temp', 
                    self.seq, self.quality, ref))

    def _cpg_methylation_average(self):
        x = []
        y = []
        window = 100
        cpg_met_scores = self.guppy_mbt[:,3]
        for n, item in enumerate([cpg_met_scores[i:i+window] for i in range(0, len(cpg_met_scores), window)]):
            x.append(n * window)
            y.append(len([i for i in item if i>180]) / window * 70000)
        return x, y

    def plot_structure(self, ref, savedir):
        sam_sum = utilities.split_mapping_and_sam_analysis(self.split_length, self.read_id, self.seq, self.quality, ref)
        lc = utilities.plot_read_structure(self.read_id, self.split_length, sam_sum)
        fig = plt.figure()
        plt.subplots_adjust(left=0.2)
        ax = fig.add_subplot()
        x, y = self._cpg_methylation_average()
        for n, item in enumerate(self.guppy_mbt[:,1]):
            if item > 120:
                ax.bar(n, item * 10000/255, width=100, color='red', zorder=3)
        ax.bar(x, y, width=100, color = 'mediumblue', zorder=0)
        ax.add_collection(lc)
        ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
        ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
        ax.autoscale()
        ax.set_ylim([-12000, 46000])
        if savedir[-1] == '/':
            savedir = savedir[:-1]
        plt.savefig(savedir + '/' + self.read_id + '.png', dpi=300)
        plt.close()

    def get_coding_coordinates(self):
        coding_start = 0
        split_length = 300
        coding_end = 13000
        positions = self.sam_summary[:,2].astype('int32')
        in_codings = [n if 0 < n < 13000 else 0 for n in positions]
        flag = 0
        start = 0 
        end = 0
        coding_starts_ends = []
        for n, pos in enumerate(in_codings):
            if flag == 0 and pos != 0:
                flag = 1
                start = n
            if flag == 1 and pos == 0:
                flag = 0
                end = n
                if end - start > 5:
                    coding_starts_ends.append((start * split_length, end * split_length))
        if flag == 1:
            end = n
            if end - start > 5:
                coding_starts_ends.append((start * split_length, end * split_length))
        return coding_starts_ends

    def get_coding_met_stats(self):
        coding_pos = self.get_coding_coordinates()
        cpg_stats = self.guppy_mbt[:,3]
        for start, end in coding_pos:
            print(np.mean(cpg_stats[start:end]))

if __name__ == '__main__':
    ref = 'rDNA_index/humRibosomal.fa'
    fast5_dir = '/mnt/data2/hpgp/HG02717_1_rDNA_fast5s_bc/workspace'
    fast5files = glob.glob(fast5_dir + '/*.fast5', recursive=True)
    reads = []
    s = 0
    e = 500
    for n, f in enumerate(fast5files):
        if n < 5:
            continue
        r = Fast5read(f, ref)
        plt.plot([*range(s,e)], r.raw[s:e])
        bps = []
        print(r.guppy_move[:50])
        for n in range(10000):
            if r.guppy_move[n]:
                bps.append(r.raw_start + n*r.raw_step)
        for i in bps:
            if s < i < e:
                plt.vlines(x=i, ymin=400, ymax=700, color='black', linewidth=0.5)
        plt.show()
        quit()
    """
        for a, b in zip(r.guppy_move, r.guppy_trace):
            print(a)
            print(b)
        quit()
    """
    if os.path.exists(savedir):
        shutil.rmtree(savedir)
    os.mkdir(savedir)
    for fast5 in fast5files:
        r = Read(fast5)
        r.plot_structure(ref, savedir)

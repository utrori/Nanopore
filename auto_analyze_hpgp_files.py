import subprocess
import shutil
import os
import glob
import itertools
import main_analysis_fastq
import fast5class
import copy
import datetime
from concurrent import futures as confu


ref = 'rDNA_index/humRibosomal.fa'


def get_basename(url_str):
    basename = url_str.split('/')[-1].split('.')[0]
    if 'R941' in basename:
        splitted = basename.split('_')
        basename = '_'.join([splitted[-2], splitted[-1]])
    return basename


def compare_time(cand_list):
    current = ''
    time2fastq = {}
    for cand in cand_list:
        dt_str = ' '.join(cand.split()[:2])
        dt_formatted = datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
        time2fastq[dt_formatted] = cand
        if not current:
            current = dt_formatted
        else:
            if dt_formatted > current:
                current = dt_formatted
    return time2fastq[current]


def analyze_fastq(filename, savename):
    ret_str = ''
    total_length = 0
    read_n = 0
    with open(filename) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            read = main_analysis_fastq.Fastqread(items, ref)
            total_length += read.length
            if read._is_this_healthy_rDNA():
                ret_str += ''.join(items)
    with open(savename, 'w') as fw:
        fw.write(ret_str)
    with open('hpgp_infos.txt', 'a') as fw:
        fw.write(filename + '\t' + str(total_length) + '\t' + str(read_n) + '\n')


def download_and_tar_fast5(fast5_download):
    subprocess.run(fast5_download, shell=True)
    subprocess.run('tar -xf *.tar', shell=True)


def download_files():
    with open('hpgp.txt') as f:
        hpgp_file = f.readlines()
    hpgp_file2 = copy.copy(hpgp_file)
    analyzed_sample_names = []
    for line in hpgp_file:
        if 'nanopore' in line and 'fast5' in line and ('TiB' in line or 'GiB' in line) and 'submissions' not in line and 'circulomics' not in line and 'circ' not in line:
            row = line.strip().split()
            if 'TiB' in line:
                continue
            if 'GiB' in line:
                if float(row[2]) > 60:
                    continue
            print(line)
            fast5_download = 'aws s3 cp s3://human-pangenomics/' + row[4] + ' /mnt/data/hpgp-data/temp'
            basename = get_basename(row[4])
            print(basename)
            sample_name = basename.split('_')[0]
            if sample_name in analyzed_sample_names:
                continue
            candidate_fastqs = []
            path_and_base = '/'.join(line.split()[4].split('/')[:-1] + [basename])
            for line2 in hpgp_file2:
                if 'fastq' in line2 and 'nanopore' in line2 and 'fast5' not in line2:
                    if path_and_base + '.' in line2 or path_and_base + '_' in line2:
                        candidate_fastqs.append(line2.strip())
            if candidate_fastqs == []:
                with open('fast5s_wo_fastq.txt', 'a') as fw:
                    fw.write(basename + '\n')
            else:
                analyzed_sample_names.append(sample_name)
                """
                fastq_downlaod = 'aws s3 cp s3://human-pangenomics/' + compare_time(candidate_fastqs).split()[4] + ' /mnt/data/hpgp-data/temp'
                subprocess.run(fastq_downlaod, shell=True)
                """
                fg = glob.glob('/mnt/data/hpgp-data/temp/*.fastq.gz')[0]
                subprocess.run('gunzip -d ' + fg, shell=True)
                rDNA_fastq_filename = '/mnt/data/hpgp-data/' + basename + '_rDNA.fastq'
                with confu.ThreadPoolExecutor(max_workers=2) as executor:
                    futures = [analyze_fastq(fg[:-3], rDNA_fastq_filename), 
                               download_and_tar_fast5(fast5_download)]
                    confu.wait(futures)
                os.remove(fg[:-3])
                dirs = glob.glob('/mnt/data/hpgp-data/temp/**/')
                out_dir = '/mnt/data/hpgp-data/' + basename + '_rDNA_fast5s'
                os.mkdir(out_dir)
                find_rDNA_fast5(rDNA_fastq_filename, dirs[0], out_dir)
                shutil.rmtree(dirs[0])


def find_rDNA_fast5(rDNA_fastq, fast5_dir, out_dir):
    if os.path.exists('/mnt/data/hpgp-data/temp_fast5s'):
        shutil.rmtree('/mnt/data/hpgp-data/temp_fast5s')
    read_ids = []
    with open(rDNA_fastq) as f:
        for items in itertools.zip_longest(*[iter(f)]*4):
            read_ids.append(items[0].split()[0][1:])
    fast5_path = glob.glob(fast5_dir + '**/*.fast5', recursive=True)
    for fast5 in fast5_path:
        subprocess.run('multi_to_single_fast5 -t 6 --recursive -i ' + fast5 + ' -s /mnt/data/hpgp-data/temp_fast5s -t 6', shell=True) 
        for split_fast5 in glob.glob('/mnt/data/hpgp-data/temp_fast5s/**/*.fast5', recursive=True):
            if fast5.split('/')[-1].split('.')[0] in read_ids:
                shutil.copy(fast5, out_dir)
        shutil.rmtree('/mnt/data/hpgp-data/temp_fast5s')


def basecall(in_dir):
    if in_dir[-1] == '/':
        in_dir = in_dir[:-1]
    subprocess.run('~/Softwares/ont-guppy_4.2.2_linux64/ont-guppy/bin/guppy_basecaller -i {} -s {}_bc -c dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac_prom.cfg --device cuda:0 --fast5_out'.format(in_dir, in_dir), shell=True)
    shutil.rmtree(in_dir)


if __name__ == '__main__':
    download_files()
    quit()
    find_rDNA_fast5('/mnt/data2/hpgp/HG02717_1_rDNA.fastq', '/mnt/data2/hpgp/08_11_20_R941_HG02717/HG02717_1/', 'HG02717_1_rDNA_fast5s')

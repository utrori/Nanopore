import os
import pandas as pd
import itertools
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
from visualize_nanopore_read import plot_read_structure


FNULL = open(os.devnull, 'w')


def align_pair(seq1, seq2):
    """Perform EMBOSS Needle for a pair of seq

    Args:
        seq1 (str): nucleotide seq
        seq2 (str): nucleotide seq
    Returns:
        score (int): alignment score
    """
    length = 4000
    with open('temp_files/temp_seq1.fa', 'w') as fw:
        fw.write('>temp1\n' + seq1[:length])
    with open('temp_files/temp_seq2.fa', 'w') as fw:
        fw.write('>temp2\n' + seq2[:length])
    subprocess.run('needle -asequence temp_files/temp_seq1.fa -bsequence '
                   'temp_files/temp_seq2.fa -outfile out.needle '
                   '-gapopen 10 -gapextend 0.5',
                   shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    with open('out.needle') as f:
        res = f.readlines()
    score = float(res[26][9:].strip())
    return score


def group_reads(seq_list):
    """Perform alignment for all the pairs in the list.

    This function is based on DP

    Args:
        seq_list (list): list of sequences. each seq is (index, seq).
    Returns:
        groups (list): a list of similar reads
    """
    groups = []
    removed = []
    for n in range(len(seq_list)):
        print(n)
        if n not in removed:
            # group is a list of sequences.
            # Each sequence is expressed by its number.
            group = [n]
            for m in range(len(seq_list)):
                if m not in removed + [n]:
                    score = align_pair(seq_list[n][1], seq_list[m][1])
                    if score > 10000:
                        group.append(m)
                        removed.append(m)
            groups.append(group)
        removed.append(n)
    return groups


if __name__ == '__main__':
    group1 = pd.read_pickle('group1.pkl')
    group2 = pd.read_pickle('group2.pkl')
    num2id1 = {}
    with open('boundary_seq1.fa') as f:
        for n, line in enumerate(f):
            if n % 3 == 0:
                num2id1[n // 3] = line.strip()[1:]
    num2id2 = {}
    with open('boundary_seq2.fa') as f:
        for n, line in enumerate(f):
            if n % 3 == 0:
                num2id2[n // 3] = line.strip()[1:]
    for n, i in enumerate(group1):
        if len(i) > 2:
            for item in i:
                plot_read_structure(num2id1[item], 200, savename='boundary1/'
                                   + str(n) + '/' + num2id1[item])
    for n, i in enumerate(group2):
        if len(i) > 2:
            try:
                os.mkdir('boundary2/' + str(n))
            finally:
                for item in i:
                    plot_read_structure(num2id2[item], 200, savename='boundary2/'
                                       + str(n) + '/' + num2id2[item])
    quit()
    seq_list1 = []
    with open('boundary_seq1.fa') as f:
        for fa in itertools.zip_longest(*[iter(f)]*3):
            seq_list1.append((fa[0].strip(), fa[1].strip()))

    seq_list2 = []
    with open('boundary_seq2.fa') as f:
        for fa in itertools.zip_longest(*[iter(f)]*3):
            seq_list2.append((fa[0].strip(), fa[1].strip()))

    groups = group_reads(seq_list1)
    print(groups)
    pd.to_pickle(groups, 'group1.pkl')
    groups = group_reads(seq_list2)
    print(groups)
    pd.to_pickle(groups, 'group2.pkl')

    quit()
    scores = []
    l = 5
    for n in range(0, len(seq_list1)):
        if n != l:
            score = align_pair(seq_list1[l][1], seq_list1[n][1])
            scores.append(score)

    sns.set()
    sns.set_style('whitegrid')
    sns.set_palette('gray')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(scores)
    plt.show()

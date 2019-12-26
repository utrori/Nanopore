from rDNA_structure_of_long_reads import easy_flag, analyze_fastq_by_header
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.collections as mc
import matplotlib.cm as cm


def plot_read_structure(header, split_length, savename=None, title=None):
    """Plot read based on split mapped sam file.

    Args:
        header (str): header name
        split_length (int): split_length of the read when it was mapped to rDNA
        savename (str): if this is specifed, the plot is saved. otherwise,
        it shows the plot.

    Returns:
        None
    """
    analyze_fastq_by_header('rDNA_reads.fastq', header, split_length)
    samfile = 'temp_files/single_split_mapped.sam'

    plot = []
    offset = 0
    rD_size = 42999
    TR_size = 13314
    with open(samfile) as f:
        samdata = f.readlines()[2:]
        for line in samdata:
            row = line.split()
            flag = int(row[1])
            # temporarily disregard multiply mapped reads
            if easy_flag(flag, 4) == 1:
                plot.append((-10000, '*'))
            else:
                if easy_flag(flag, 16) != 1:
                    plot.append((int(row[3]) + offset, '+'))
                else:
                    plot.append((int(row[3]) + offset, '-'))
    read_num = len(plot)

    x = np.linspace(0, split_length * (read_num - 1), read_num)
    r_coord = [i[0] for i in plot]
    vertical_lines = []
    for n in range(read_num):
        if plot[n][1] == '+':
            line_element = [(x[n], r_coord[n]),
                            (x[n] + split_length, r_coord[n] + split_length)]
        elif plot[n][1] == '-':
            line_element = [(x[n], r_coord[n]),
                            (x[n] + split_length, r_coord[n] - split_length)]
        else:
            line_element = [(x[n], r_coord[n]),
                            (x[n] + split_length, r_coord[n])]
        vertical_lines.append(line_element)
    vertical_lines.append([[0, offset], [split_length * read_num, offset]])
    vertical_lines.append([[0, offset + TR_size],
                          [split_length * read_num, offset + TR_size]])
    vertical_lines.append([[0, offset + rD_size],
                          [split_length * read_num, offset + rD_size]])

    lw = [1 for i in range(read_num)]
    lw.extend([0.5 for i in range(3)])
    ls = ['solid' for i in range(read_num)]
    ls.extend(['dashed' for i in range(3)])
    cl = []
    for coord in r_coord:
        if coord == -10000:
            cl.append('black')
        else:
            cl.append(cm.hsv((coord - offset) / rD_size))
    cl.extend(['black' for i in range(3)])
    lc = mc.LineCollection(vertical_lines, linewidths=lw,
                           linestyles=ls, colors=cl)

    fig = plt.figure()
    plt.subplots_adjust(left=0.2)
    ax = fig.add_subplot()
    ax.add_collection(lc)
    ax.autoscale()
    ax.set_yticks((-10000, 0, 10000, 20000, 30000, 40000))
    ax.set_yticklabels(('unmapped', 0, 10000, 20000, 30000, 40000))
    if title:
        ax.set_title(title)
    if not savename:
        plt.show()
    else:
        plt.savefig(savename, dpi=300)
    plt.close()


if __name__ == '__main__':
    split_length = 200
    #header = sys.argv[1]
    headers = []
    header = '@1e384cf6-b197-4642-88a9-7f120013b16d'
    plot_read_structure(header, split_length)
    quit()
    with open('reads_visualized.txt') as f:
        for n, line in enumerate(f):
            if n % 2 == 0:
                header = line.strip().split()[1]
                headers.append(header)
                savename = 'reads_visualized/figure_' + str(n // 2) + '.png'
                plot_read_structure(header, split_length, savename)
            else:
                continue

def print_genbank(seq):
    """Print a sequence as a formatted genbank format

    Args:
        seq (str): sequence string
    Returns:
        gen (str): genbank formatted sequence
    """
    gen = ''
    for n in range((len(seq) - 1) // 60 + 1):
        pos = n * 60 + 1
        digit = len(str(pos))
        gen += ' ' * (9 - digit) + str(pos)
        for m in range(6):
            gen += ' ' + seq[pos+m*10-1:pos+m*10+9]
        gen += '\n'
    return gen

if __name__ == '__main__':
    header = '58ef3b9-bbd2-4990-97ae-d82085669e2d'
    with open('rDNA_reads.fastq') as f:
        flag = 0
        for line in f:
            if flag == 0:
                if header in line:
                    flag = 1
            else:
                print(print_genbank(line))
                quit()

with open('chr21.sam') as f:
    f.readline()
    f.readline()
    for line in f:
        if line.split()[2] != '*':
            print(line.split('_')[1].split()[0])

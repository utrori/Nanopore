def read_genbank(filename):
    features = []
    with open(filename) as f:
        for line in f:
            if 'misc_feature' in line:
                row = line.split()
                position = row[1].split('..')
            if '/label=' in line:
                label = line.split('\"')[1].split('\"')[0]
                features.append((label, position))
        return features


if __name__ == '__main__':
    filename = '/home/yutaro/Downloads/Human_rDNA_consensus.ape'
    features = read_genbank(filename)
    print features

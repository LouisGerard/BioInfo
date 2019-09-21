def probability(kmer, profile):
    return sum([profile[i][nuc] for i, nuc in enumerate(kmer)])

def most_probable_kmer(text, profile):
    k = len(profile)
    best_kmer = ''
    best_p = 0
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = probability(kmer, profile)
        if p > best_p:
            best_p = p
            best_kmer = kmer
    return best_kmer

if __name__ == '__main__':
    fname = 'test.txt'
    with open(fname) as f:
        gene = f.readline().rstrip('\n')
        k = int(f.readline().rstrip('\n'))
        profile = [{} for i in range(k)]
        for nuc in ['A', 'C', 'G', 'T']:
            line = f.readline().rstrip('\n')
            for i, p in enumerate(line.split()):
                p = float(p)
                profile[i][nuc] = p

    print(most_probable_kmer(gene, profile))

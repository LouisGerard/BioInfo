import itertools
from collections import Counter
from math import log

from week2 import neighbors, hamming


def motif_in_text(text, motif, d):
    k = len(motif)
    for i in range(len(text) - k + 1):
        if hamming(text[i:i+k], motif) <= d:
            return True
    return False


def motif_in_all_texts(texts, motif, d):
    for text in texts:
        if not motif_in_text(text, motif, d):
            return False
    return True


def motif_enum(texts, k, d):
    results = set()
    full_text = ''.join(texts)
    for i in range(len(full_text) - k + 1):
        kmer = full_text[i:i+k]
        for d_ in range(d+1):
            for n in neighbors(kmer, d_):
                if motif_in_all_texts(texts, n, d):
                    results.add(n)
    return results


def profile(matrix):
    cols = len(matrix[0])
    counts = [Counter() for _ in range(cols)]
    for row in matrix:
        for i, nuc in enumerate(row):
            counts[i][nuc] += 1
    sums = [sum(counts[i].values()) for i in range(cols)]
    return [{k: v / sums[i] for k, v in counts[i].items()} for i in range(cols)]


def entropy(matrix):
    p = profile(matrix)
    result = 0
    for col in p:
        for v in col.values():
            result -= log(v, 2) * v
    return result


def motifs(texts, pattern):
    k = len(pattern)
    result = []
    dists = []
    for t in texts:
        best_kmer = ''
        min_dist = float('inf')
        for i in range(len(t) - k + 1):
            kmer = t[i:i+k]
            dist = hamming(pattern, kmer)
            if dist < min_dist:
                min_dist = dist
                best_kmer = kmer
        result.append(best_kmer)
        dists.append(min_dist)
    return result, dists


def brute_force_median_string(texts, k):
    best_kmer = ''
    min_dist = float('inf')
    for i in itertools.product('ACTG', repeat=k):
        kmer = ''.join(i)
        dist = sum(motifs(texts, kmer)[1])
        if dist < min_dist:
            min_dist = dist
            best_kmer = kmer
    return best_kmer


if __name__ == '__main__':
    fname = '/home/louis/Documents/Bio_courses/dataset_158_9.txt'
    with open(fname) as f:
        params = f.readline().rstrip('\n')
        k = int(params)
        genes = [i.rstrip('\n') for i in f.readlines() if len(i) > 1]

    print(brute_force_median_string(genes, k))

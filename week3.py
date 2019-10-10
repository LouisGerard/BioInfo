import itertools
from collections import Counter
from math import log
from functools import reduce
import operator

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
    sums = [sum([counts[i][j] for j in 'ACTG']) for i in range(cols)]
    return [{k: counts[i][k] / sums[i] for k in 'ACTG'} for i in range(cols)]


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
    return best_kmer, min_dist


def probability(kmer, profile):
    return reduce(operator.mul, [profile[i][nuc] for i, nuc in enumerate(kmer)])


def most_probable_kmer(text, profile):
    k = len(profile)
    best_kmer = text[:k]
    best_p = probability(best_kmer, profile)
    for i in range(1, len(text) - k + 1):
        kmer = text[i:i+k]
        p = probability(kmer, profile)
        if p > best_p:
            best_p = p
            best_kmer = kmer
    return best_kmer


def consensus(profile):
    result = ''
    for i in profile:
        result += max(i, key=i.get)
    return result


def score(profile):  # sum profile is equivalent to sum hamming
    result = 0
    c = consensus(profile)
    for i, c_ in zip(profile, c):
        result += sum([v for k, v in i.items() if k != c_])
    return result


def greedy_motif_search(texts, k, start_from_1=False):
    missing_fn = Counter.__missing__
    if start_from_1:
        Counter.__missing__ = lambda self, x: 1

    t = len(texts)

    best_motifs = [text[:k] for text in texts]
    best_score = float('inf')
    for i in range(len(texts[0]) - k + 1):
        motifs = [texts[0][i:i+k]]
        for j in range(1, t):
            p = profile(motifs)
            motifs.append(most_probable_kmer(texts[j], p))
        p = profile(motifs)
        s = score(p)

        if s < best_score:
            best_score = s
            best_motifs = motifs

    Counter.__missing__ = missing_fn
    return best_motifs


if __name__ == '__main__':
    texts = [
        'CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
        'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
        'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG',
    ]
    print('\n'.join(brute_force_median_string(texts, 7)))

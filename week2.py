import itertools
from collections import Counter

from matplotlib import pyplot as plt

from week1 import complement


def skew(text):
    result = [0]
    for i in range(len(text)):
        if text[i] == 'C':
            result.append(result[i] - 1)
        elif text[i] == 'G':
            result.append(result[i] + 1)
        else:
            result.append(result[i])
    return result


def min_skew(text):
    s = skew(text)
    min_ = min(s)
    return [i for i in range(len(s)) if s[i] == min_]


def hamming(text1, text2):
    cnt = 0
    for i1, i2 in zip(text1, text2):
        if i1 != i2:
            cnt += 1
    return cnt


def approximate_pattern_matching(text, needle, d):
    k = len(needle)
    results = []
    for i in range(len(text) - k + 1):
        if hamming(needle, text[i:i + k]) <= d:
            results.append(i)
    return results


def count(text, needle, d):
    k = len(needle)
    results = 0
    for i in range(len(text) - k + 1):
        if hamming(needle, text[i:i + k]) <= d:
            results += 1
    return results


def neighbors(pattern, d):
    nucleotides = {'A', 'C', 'T', 'G'}
    k = len(pattern)

    for pos in itertools.combinations(range(k), d):
        neighbors_nuc = [nucleotides - {pattern[j]} for j in pos]
        for variation in itertools.product(*neighbors_nuc):
            variation_i = 0
            neighbor = ''
            for pattern_i in range(k):
                if pattern_i in pos:
                    neighbor += variation[variation_i]
                    variation_i += 1
                else:
                    neighbor += pattern[pattern_i]
            yield neighbor


def most_frequents(text, k, d, count_complement=True):
    counts = Counter()

    for i in range(len(text) - k + 1):
        kmer = text[i:i + k]
        counts[kmer] += 1
        if count_complement:
            counts[complement(kmer)] += 1

    added = Counter()
    for kmer, count in counts.items():
        for d_ in range(d):
            for neighbor in neighbors(kmer, d+1):
                added[neighbor] += count
    counts += added

    max_ = max(counts.values())
    return [kmer for kmer, count in counts.items() if count == max_]


def visu_skew(text):
    s = skew(text)
    plt.plot(s)
    plt.show()


if __name__ == '__main__':
    print(len(list(neighbors('TGCAT', 1))))

import random
import sys
from collections import Counter

from week3 import profile, most_probable_kmer, score


def randomized_motif_search(texts, k):
    missing_fn = Counter.__missing__
    Counter.__missing__ = lambda self, x: 1

    init = [random.randint(0, len(t) - k) for t in texts]
    motifs = [t[init[i]:init[i]+k] for i, t in enumerate(texts)]
    best_score = float('inf')
    while True:
        p = profile(motifs)
        s = score(p)
        if s < best_score:
            best_score = s
            motifs = [most_probable_kmer(text, p) for text in texts]
        else:
            Counter.__missing__ = missing_fn
            return motifs


def sampled_randomized_motif_search(texts, k, n=1000):
    result = None
    best_score = float('inf')
    for _ in range(n):
        motifs = randomized_motif_search(texts, k)
        s = score(profile(motifs))
        if s <= best_score:
            best_score = s
            result = motifs
    return result


if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        params = f.readline().rstrip('\n')
        k = int(params.split(maxsplit=1)[0])
        texts = [line.rstrip('\n') for line in f]

    result = sampled_randomized_motif_search(texts, k)

    print('\n'.join(result))

import random
import sys
from collections import Counter

from tqdm import tqdm

from week3 import profile, most_probable_kmer, score, probability


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


def sampled_randomized_motif_search(*args, n=1000, search_fn=randomized_motif_search):
    result = None
    best_score = float('inf')
    for _ in tqdm(range(n)):
        motifs = search_fn(*args)
        s = score(profile(motifs))
        if s <= best_score:
            best_score = s
            result = motifs
    print(best_score)
    return result


def profile_random_generation(text, p):
    kmers = []
    probas = []

    k = len(p)

    for i in range(len(text) - k + 1):
        kmers.append(text[i:i+k])
        probas.append(probability(kmers[-1], p))

    return random.choices(kmers, probas)[0]


def gibbs_motif_search(texts, k, n=100):
    missing_fn = Counter.__missing__
    Counter.__missing__ = lambda self, x: 1

    t = len(texts)
    init = [random.randint(0, len(t) - k) for t in texts]
    motifs = [t[init[i]:init[i]+k] for i, t in enumerate(texts)]

    best_score = float('inf')
    best_motifs = None

    for _ in range(n):
        i = random.randint(0, t-1)
        p = profile(motifs[:i] + motifs[i+1:])
        motifs[i] = profile_random_generation(texts[i], p)
        s = score(p)
        if s < best_score:
            best_score = s
            best_motifs = motifs[:]

    Counter.__missing__ = missing_fn
    return best_motifs


def iterative_xor(p):
    xor = 2 * p - 2 * p * p
    nothing = 1 - p
    for _ in range(8):
        nothing *= (1 - p)
        xor += (nothing - xor) * p
    nothing *= (1 - p)
    return 1 - nothing - xor


if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        params = f.readline().rstrip('\n').split()
        k = int(params[0])
        n = int(params[2])
        texts = [line.rstrip('\n') for line in f]

    print('\n'.join(sampled_randomized_motif_search(texts, k, n, n=100, search_fn=gibbs_motif_search)))

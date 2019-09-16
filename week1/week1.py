from collections import Counter

def pattern_count(text, needle):
    count = 0
    k = len(needle)
    for i in range(len(text) - k + 1):
        if text[i:i+k] == needle:
            count += 1
    return count

def all_patterns_counts(text, k):
    counts = []
    for i in range(len(text) - k + 1):
        counts.append(pattern_count(text, text[i:i+k]))
    return counts

def frequent_words(text, k):
    counts = all_patterns_counts(text, k)
    max_ = max(counts)

    result = set()
    for i, c in enumerate(counts):
        if c == max_:
            result.add(text[i:i+k])
    return result

def complement(text):
    translation = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
    }
    result = ''
    for i in reversed(text):
        result += translation[i]
    return result

def substr_pos(text, needle):
    result = []
    k = len(needle)
    for i in range(len(text) - k + 1):
        if text[i:i+k] == needle:
            result.append(i)
    return result

def find_clumps(text, k, L=500, t=3):
    pos = {}
    result = set()
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]

        # print('-'*50)
        # print(kmer)
        if kmer not in pos:
            # print('not in pos')
            pos[kmer] = [i]
            continue
        if kmer in result:
            # print('in result')
            continue

        to_del = []
        for j, p in enumerate(pos[kmer]):
            if i - p + k > L:
                to_del.append(j)
                continue
        # print('to del :', to_del)
        for j in reversed(to_del):
            pos[kmer].pop(j)
        if len(pos[kmer]) + 1 >= t:
            # print('add result')
            result.add(kmer)
            pos.pop(kmer)
            continue
        # print('add pos')
        pos[kmer].append(i)
        # print('pos :', pos[kmer])
    return result

if __name__ == '__main__':
    # with open('E_coli.txt') as f:
    #     gene = f.readline().rstrip('\n')
    # k = 9
    # L = 500
    # t = 3

    print(substr_pos('AAACATAGGATCAAC', 'AA'))
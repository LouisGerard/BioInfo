def skew(text):
    result = [0]
    for i in range(len(text)):
        if text[i] == 'C':
            result.append(result[i]-1)
        elif text[i] == 'G':
            result.append(result[i]+1)
        else:
            result.append(result[i])
    return result


def min_skew(text):
    s = skew(text)
    min_ = min(s)
    return [i for i in range(len(s)) if s[i] == min_]


if __name__ == '__main__':
    with open('/home/louis/Documents/Bio_courses/dataset_7_6.txt') as f:
        gene = f.readline().rstrip('\n')
    print(' '.join(map(str, min_skew(gene))))

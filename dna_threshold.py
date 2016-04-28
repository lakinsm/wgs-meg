#!/usr/bin/env python3

import sys

dna_rep = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
ambig = ('N', 'Y', 'W', 'K', 'I', 'B', 'D')

phredmap = {
    '!': 0, '\"': 1, '#': 2, '$': 3, '%': 4, '&': 5, '\'': 6,
    '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13, '/': 14,
    '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20, '6': 21, '7': 22,
    '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30,
    '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38,
    'H': 39, 'I': 40, 'J': 41, 'K': 42, 'L': 43, 'M': 44, 'N': 45, 'O': 46,
    'P': 47, 'Q': 48, 'R': 49, 'S': 50, 'T': 51, 'U': 52, 'V': 53, 'W': 54,
    'X': 55, 'Y': 56, 'Z': 57, '[': 58, '\\': 59, ']': 60, '_': 61, '`': 62,
    'a': 63, 'b': 64, 'c': 65, 'd': 66, 'e': 67, 'f': 68, 'g': 69, 'h': 70,
    'i': 71, 'j': 72, 'k': 73, 'l': 74, 'm': 75, 'n': 76, 'o': 77, 'p': 78,
    'q': 79, 'r': 80, 's': 81, 't': 82, 'u': 83, 'v': 84, 'w': 85, 'x': 86,
    'y': 87, 'z': 88, '{': 89, '|': 90, '}': 91, '~': 92
}

thresh = int(sys.argv[2])


def parse_fastq(infile):
    while True:
        header = infile.readline()[1:].rstrip()
        seq = infile.readline().strip()
        if not header or not seq:
            break
        _ = infile.readline()
        scores = infile.readline()
        yield header, seq, scores

if __name__=='__main__':
    with open(sys.argv[1], 'r') as inf:
        for h, seq, score in parse_fastq(inf):
            outseq = [x if phredmap[score[e]] >= thresh-1 else 'N' for e, x in enumerate(seq)]
            sys.stdout.write('>'+h+'\n'+''.join(outseq)+'\n')












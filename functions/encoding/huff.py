from heapq import heappush, heappop, heapify
from collections import Counter
import os
import numpy as np
import scipy.io

def encode(symb2freq):
    """Huffman encode the given dict mapping symbols to weights"""
    heap = [[wt, [sym, ""]] for sym, wt in symb2freq.items()]
    heapify(heap)
    while len(heap) > 1:
        lo = heappop(heap)
        hi = heappop(heap)
        for pair in lo[1:]:
            pair[1] = '0' + pair[1]
        for pair in hi[1:]:
            pair[1] = '1' + pair[1]
        heappush(heap, [lo[0] + hi[0]] + lo[1:] + hi[1:])
    return sorted(heappop(heap)[1:], key=lambda p: (len(p[-1]), p))

def main():

    fname = os.sys.argv[1]

    mat = scipy.io.loadmat(fname)
    mat_data = mat['q_data'].tolist()
    data = []
    for e in mat_data[0]:
        # print(e)
        data.append(int(e))

    huffify(data)
    return

def test():
    return 42

def huffify(data):

    symb2freq = Counter(data)

    huff = encode(symb2freq)
    # print ("Symbol\tWeight\tHuffman Code")
    num_bits = 0
    huff_overhead = 0
    for p in huff:
        # print ("%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1]))

        num_bits += symb2freq[p[0]]* len(p[1])
        huff_overhead += 2*16

    # print((huff_overhead+num_bits)/(len(data)*16))
    # print(huff_overhead+num_bits)
    return num_bits

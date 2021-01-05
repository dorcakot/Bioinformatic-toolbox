import numpy
from Bio import pairwise2

def create_matrix(seq1, seq2):
    """
    Using dynamic programming, function creates matrix containing edit distance values for each i-long and j-long
    refixes of seq1 and seq2 at positions [i, j].
    :param seq1, seq2: sequences for evaluation
    :return: matrix D with partial edit distance values
    """
    D = numpy.zeros((len(seq1) + 1, len(seq2) + 1), dtype=numpy.int32)
    for i in range(len(seq1) + 1):
        D[i, 0] = i
    for j in range(len(seq2) + 1):
        D[0, j] = j
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if seq1[i - 1] == seq2[j - 1]:
                D[i, j] = D[i - 1, j - 1]
            elif seq1[i - 1] != seq2[j - 1]:
                D[i, j] = min(D[i - 1, j - 1] + 1, D[i - 1, j] + 1, D[i, j - 1] + 1)
    return D


def edit_distance(seq1, seq2):
    """
    Function takes two sequences, calls function create_matrix, which calculates edit distance for each i-th and j-th
    prefixes of given sequences and return edit distance of these sequences, which is found in the D[|seq1|, |seq2|] position.
    :param seq1, seq2: sequences for evalutation
    :return: edit distance of two given sequences
    """
    D = create_matrix(seq1, seq2)
    return D[len(seq1)][len(seq2)]


def all_alignments(seq1, seq2):
    """
    This function takes two seqneces and computes their all posiible optimal alignments using pairwise2 library.
    :param seq1, seq2: sequences to be evaluated
    :return: all possible optimal alignments of seq1, seq2
    """
    alignments = []
    for alignment in pairwise2.align.globalms(seq1, seq2, 0, -1, -1, -1):
        alignments.append(alignment[0] +"\n"+ alignment[1])
    return alignments

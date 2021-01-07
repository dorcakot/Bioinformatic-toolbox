from fasta import Fasta
from typing import Union


def hd(seq1: Union[str, Fasta], seq2: Union[str, Fasta]) -> int:
    """
    Function calculates Hamming distance of two sequences of the same length. When sequences are not of the same length,
    an exception is raised.
    :param seq1, seq2: sequences to be evaluated
    :return: Hamming distance value for two given sequences
    """
    if len(seq1) != len(seq2):
        raise Exception('Sequences have different length!')
    ham_dist=0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            ham_dist += 1
    return ham_dist

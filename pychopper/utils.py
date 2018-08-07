
import numpy as np


BASES = ('A', 'T', 'G', 'C')


def seq_to_array(seq):
    """Convert sequence to numpy array of integers.

    :param seq: DNA sequence.
    :returns: Sequence as an array of integers.
    :rtype: numpy.array
    """
    result = np.zeros(len(seq), dtype='i4')
    for i, base in enumerate(seq):
        result[i] = BASES.index(base)
    return result


def random_seq(length):
    """ Generate random sequence of specified length.
    :param length: Length of seqeucne to generate.
    :returns: Random DNA sequence,
    :rtype: str
    """
    return ''.join([BASES[b] for b in np.random.random_integers(0, len(BASES) - 1, size=length)])

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import parasail

BASES = ('A', 'T', 'G', 'C')

DEFAULT_ALIGN_PARAMS = {'match': 1,
                        'mismatch': -1,
                        'gap_open': 2,
                        'gap_extend': 1}


def random_seq(length):
    """ Generate random sequence of specified length.
    :param length: Length of sequence to generate.
    :returns: Random DNA sequence.
    :rtype: str
    """
    return ''.join([BASES[b] for b in np.random.random_integers(0, len(BASES) - 1, size=length)])


def pair_align(reference, query, params=DEFAULT_ALIGN_PARAMS):
    """ Perform pairwise local alignment using scikit-bio.
    :param reference: Reference sequence.
    :param query: Query sequence.
    :param params: Alignment parameters in a dictionary.
    :returns: Alignments in scikit-bio format.
    :rtype: list of tuples
    """

    subs_mat = parasail.matrix_create("ACGT", params['match'], params['mismatch'])
    aln = parasail.sw_striped_32(reference, query, params['gap_open'], params['gap_extend'], subs_mat)

    return aln


def score_cutoff(barcode, aln_params, target_length, percentile, nr_samples):
    """ Calculate a score cutoff for a barcode by aligning it to random sequences.

    :param barcode: The query sequence.
    :param aln_params: Alignmment parameters.
    :param target_length: Length of simulated random sequences.
    :param percentile: The percentile of score distribution used as score cutoff.
    :param nr_samples: Number of random sequences to align against.
    :returns: The calculated score cutoff.
    :rtype: float
    """
    score_iter = (pair_align(random_seq(target_length), barcode, params=aln_params)
                  for i in range(nr_samples))
    null_scores = []
    for aln in score_iter:
        null_scores.append(aln.score)
    # cutoff = np.percentile(null_scores, percentile) + np.std(null_scores) * 2
    cutoff = np.percentile(null_scores, percentile)
    return cutoff


def process_alignment(aln):
    """ Process a Bio.pairwise2 alignment, extracting score, start and end.
    :param aln: Alignment.
    :return: A dictionary with score, ref_start and ref_end.
    :rtype: dict
    """
    res = {}
    res['score'] = aln.score
    res['ref_start'] = None
    res['ref_end'] = None
    return res


def seq_detect(reference, query, score_cutoff, params=DEFAULT_ALIGN_PARAMS):
    """ Detect a query in a reference base on score cutoff and uniqueness.
    :param reference: Reference sequence.
    :param query: Query sequence.
    :param score_cutoff: Minimum alignment score.
    :param params: Alignment parameters in a dictionary.
    """
    alns = pair_align(str(reference), query, params)

    res = process_alignment(alns)
    if res['score'] < score_cutoff:
        return None
    else:
        return res


if __name__ == '__main__':

    a = "ACTTGCCTGTCGCTCTATCTTCNNNNNNTTTTTTTTTTTTTTTTTTTTVN"
    b = "TTTTTTTTTTTTTTTTTTTT"
    c = "ACTTGCCTGTCGCTCTATCTTC"

    cut_b = score_cutoff(b, DEFAULT_ALIGN_PARAMS, 200, 95, 100)
    cut_c = score_cutoff(c, DEFAULT_ALIGN_PARAMS, 200, 95, 100)
    res_b = seq_detect(a, b, cut_b)
    res_c = seq_detect(a, c, cut_c)
    print(res_b)
    print(res_c)

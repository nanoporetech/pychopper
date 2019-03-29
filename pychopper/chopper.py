# -*- coding: utf-8 -*-

from collections import OrderedDict
import numpy as np
from pychopper import seq_utils as seu
from pychopper import seq_detect


def load_barcodes(in_fas):
    """ Load barcodes from fasta and parse names into info structure.

    :param in_fas: Input fasta file.
    :return: Ordered dictionary of barcodes.
    :rtype: OrderedDict
    """
    res = OrderedDict()
    for sr in seu.read_seq_records(in_fas):
        _seq_record_to_barcode(sr, res)
    return res


def _seq_record_to_barcode(sr, res):
    """ Convert barcode name into barcode structure """
    name, bc_number = sr.id.split('|')
    barcode = OrderedDict([('name', name), ('number', bc_number),
                           ('seq', sr), ('score_cutoff', None)])
    if name not in res:
        res[name] = [barcode]
    else:
        res[name].append(barcode)
    return res


def calculate_score_cutoffs(barcodes, aln_params, target_length, percentile, nr_samples, heu=False):
    """ Calculate score cutoffs for all barcodes.
    :param barcode: group dictionary.
    :param target_length: Length of target sequence.
    :param percentile: Percentile corresponding to the calculated score cutoff.
    :param nr_samples: Number of samples to take.
    :param heu: Dummy mode.
    :return: The barcode groups with calculated cutoffs.
    :rtype: OrderedDict
    """
    for group, bcs in barcodes.items():
        for bc in bcs:
            if not heu:
                bc['score_cutoff'] = seq_detect.score_cutoff(
                    str(bc['seq']), aln_params, target_length, percentile, nr_samples)
            else:
                bc['score_cutoff'] = 0
    return barcodes


def score_barcode(seq, barcode, aln_params):
    """ Score a single barcode against a target sequence.
    :param seq: Target sequence.
    :param barcode: Barcode.
    :param aln_params: Alignment parameters.
    :return: True if alignment passes score cutoff.
    :rtype: Bool
    """
    aln = seq_detect.pair_align(
        seq, str(barcode['seq'].seq), params=aln_params)

    return (aln, aln.score >= barcode['score_cutoff'], aln.score)


def _heuristic_scoring(bc1_start_fwd, bc1_end_fwd, bc2_start_fwd, bc2_end_fwd, bc1_start_rev, bc1_end_rev, bc2_start_rev, bc2_end_rev, limit):
    """ Use heuristic scoring to determine read orientation. """
    if ((bc1_start_fwd - max(bc1_end_fwd, bc1_start_rev, bc1_end_rev)) / bc1_start_fwd >= limit) and ((bc2_end_fwd - max(bc2_start_fwd, bc2_start_rev, bc2_end_rev)) / bc2_end_fwd >= limit):
        return 'fwd_match'
    if ((bc1_start_rev - max(bc1_start_fwd, bc1_end_fwd, bc1_end_rev)) / bc1_start_rev >= limit) and ((bc2_end_rev - max(bc2_start_fwd, bc2_end_fwd, bc2_start_rev)) / bc2_end_rev >= limit):
        return 'rev_match'
    return None


def score_barcode_group(reference, target_length, barcode_group, barcodes, aln_params, heu=False, heu_limit=None):
    """ Score barcode group against target sequence and its reverse complement. Only unambigous matches
    are scored as good.

    :param reference: Target sequence.
    :param target_length: Length of sequence to consider at each end.
    :param barcode_group: A barcode group with two barcodes.
    :param aln_params: Alignment parameters.
    :param heu: Use heuristic scoring.
    :param heu_limit: Heuristic scoring stringency.
    :return: 'fwd_match' for forward match, 'rev_match' for reverese match and None for no match.
    :rtype: str or None.
    """
    bc1, bc2 = barcodes

    fwd_match = np.array([[True, False],
                          [False, True],
                          [False, False],
                          [False, False]], dtype=bool)

    rev_match = np.array([[False, False],
                          [False, False],
                          [True, False],
                          [False, True]], dtype=bool)

    target = str(reference)
    sm = np.zeros((4, 2), dtype=bool)
    aln_start, pass_start, bc1_start_fwd = score_barcode(
        target[:target_length], bc1, aln_params=aln_params)
    aln_end, pass_end, bc1_end_fwd = score_barcode(
        target[-target_length:], bc1, aln_params=aln_params)
    sm[0, 0], sm[0, 1] = pass_start, pass_end

    aln_start, pass_start, bc2_start_fwd = score_barcode(
        target[:target_length], bc2, aln_params=aln_params)
    aln_end, pass_end, bc2_end_fwd = score_barcode(
        target[-target_length:], bc2, aln_params=aln_params)
    sm[1, 0], sm[1, 1] = pass_start, pass_end

    target = seu.reverse_complement(target)
    aln_start, pass_start, bc1_start_rev = score_barcode(
        target[:target_length], bc1, aln_params=aln_params)
    aln_end, pass_end, bc1_end_rev = score_barcode(
        target[-target_length:], bc1, aln_params=aln_params)
    sm[2, 0], sm[2, 1] = pass_start, pass_end

    aln_start, pass_start, bc2_start_rev = score_barcode(
        target[:target_length], bc2, aln_params=aln_params)
    aln_end, pass_end, bc2_end_rev = score_barcode(
        target[-target_length:], bc2, aln_params=aln_params)
    sm[3, 0], sm[3, 1] = pass_start, pass_end

    nr_hits = np.sum(sm)
    score_stats = OrderedDict()

    bc1_name = bc1['seq'].name
    bc2_name = bc2['seq'].name

    score_stats[bc1_name + "_start_fwd"] = bc1_start_fwd
    score_stats[bc1_name + "_end_fwd"] = bc1_end_fwd

    score_stats[bc2_name + "_start_fwd"] = bc2_start_fwd
    score_stats[bc2_name + "_end_fwd"] = bc2_end_fwd

    score_stats[bc1_name + "_start_rev"] = bc1_start_rev
    score_stats[bc1_name + "_end_rev"] = bc1_end_rev

    score_stats[bc2_name + "_start_rev"] = bc2_start_rev
    score_stats[bc2_name + "_end_rev"] = bc2_end_rev

    if heu:
        match = _heuristic_scoring(bc1_start_fwd, bc1_end_fwd, bc2_start_fwd,
                                   bc2_end_fwd, bc1_start_rev, bc1_end_rev, bc2_start_rev, bc2_end_rev, heu_limit)
        return match, nr_hits, score_stats

    if np.all(sm == fwd_match):
        return 'fwd_match', nr_hits, score_stats
    if np.all(sm == rev_match):
        return 'rev_match', nr_hits, score_stats
    return None, nr_hits, score_stats


def score_barcode_groups(reference, barcodes, target_length, aln_params, heu=False, heu_limit=None):
    """ Score multiple barcode groups.

    :param reference: Target sequence.
    :param barcodes: A barcode groups.
    :param target_length: Length of sequence to consider at each end.
    :param aln_params: Alignment parameters.
    :param heu: Use heuristic scoring.
    :param heu_limit: Heuristic scoring stringency.
    :return: An ordered dictionary of results.
    :rtype: OrderedDict.
    """
    res = OrderedDict()
    for bc_group, barcodes in barcodes.items():
        res[bc_group] = score_barcode_group(
            reference.seq, target_length, bc_group, barcodes, aln_params, heu, heu_limit)
    return res

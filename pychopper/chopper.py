# -*- coding: utf-8 -*-

from collections import OrderedDict
import numpy as np
import itertools
from pychopper import seq_utils as seu
from pychopper import hmmer_backend, edlib_backend
from pychopper.common_structures import Segment, Seq
from pychopper.alignment_hits import process_hits


def load_primers(in_fas):
    """ Load primers from fasta and parse names into info structure.

    :param in_fas: Input fasta file.
    :return: Ordered dictionary of barcodes.
    :rtype: OrderedDict
    """
    pass


def _build_segments(hits, config):
    segments = []
    if len(hits) == 0:
        return tuple(segments)
    for s in zip(hits, hits[1:]):
        strand = None
        seg_len = 0
        c = (s[0].Query, s[1].Query)
        if c in config:
            strand = config[c]
            seg_len = s[1].RefStart - s[0].RefEnd
        segments.append(Segment(s[0].RefEnd, s[1].RefStart, strand, seg_len))
    return tuple(segments)


def analyse_hits(hits, config):
    segments = _build_segments(hits, config)
    nr = len(segments)

    if nr == 0:
        return (), (), 0

    M = np.zeros((2, nr), dtype=int)
    B = dict()

    M[1, 0] = segments[0].Len

    for j in range(1, nr):
        for i in range(2):
            if i == 0:
                M[i, j] = M[0, j - 1]
                B[i, j] = (0, j - 1)
                if M[1, j - 1] > M[0, j - 1]:
                    M[i, j] = M[1, j - 1]
                    B[i, j] = (1, j - 1)
            elif i == 1:
                M[i, j] = M[0, j - 1] + segments[j].Len
                B[i, j] = (0, j - 1)

    tlen = np.argmax(M[:, nr - 1])
    valid_segments = []

    p = (tlen, nr - 1)
    while True:
        if p[0] == 1:
            s = segments[p[1]]
            # Filter out invalid segments which can be part of the 
            # solution:
            if s.Len > 0:
                valid_segments.append(s)
        if p[1] == 0:
            break
        p = B[p]

    return tuple(valid_segments), hits, tlen


def segments_to_reads(read, segments):
    for s in segments:
        sr_seq = read.Seq[s.Start:s.End]
        sr_name = "{}:{}|".format(s.Start, s.End) + read.Name + " strand=" + s.Strand
        sr = Seq(sr_name, sr_seq, read.Qual[s.Start:s.End] if read.Qual is not None else None)
        if s.Strand == '-':
            sr = seu.revcomp_seq(sr)
        yield sr


def chopper_phmm(read, phmm_file, config, cutoff, threads):
    hits = hmmer_backend.find_locations(read, phmm_file, E=cutoff, threads=threads)
    hits = process_hits(hits, cutoff)
    return analyse_hits(hits, config)


def chopper_edlib(read, primers, config, cutoff, threads):
    pass

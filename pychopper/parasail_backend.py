# -*- coding: utf-8 -*-

import parasail
import re
from pychopper.common_structures import Hit

# Fixme: allow setting this from the command line?
DEFAULT_ALIGN_PARAMS = {'match': 1,
                        'mismatch': -2,
                        'gap_open': 1,
                        'gap_extend': 1}
DEFAULT_SUBS_MAT = parasail.matrix_create("ACGT", DEFAULT_ALIGN_PARAMS['match'], DEFAULT_ALIGN_PARAMS['mismatch'])
re_split_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)")


def first_cigar(cigar):
    """Extract details of the first operation in a cigar string."""
    m = re.search(re_split_cigar, cigar)
    return int(m.group('len')), m.group('op')


def process_alignment(aln, query, query_name, aln_params):
    """ Process an alignment, extracting score, start and end.
    """
    res = {}
    res['query_name'] = query_name
    res['score'] = aln.score
    max_score = float(aln_params['match'] * len(query))
    res['norm_score'] = (max_score - aln.score) / max_score
    fo = first_cigar(aln.cigar.decode.decode())
    res['ref_start'] = aln.cigar.beg_ref
    res['ref_end'] = aln.end_ref + 1
    res['query_start'] = aln.cigar.beg_query
    res['query_end'] = aln.end_query + 1
    if fo[1] == 'I':
        res['query_start'] += fo[0]
    if fo[1] == 'D':
        res['ref_start'] += fo[0]
    return res


def pair_align(reference, query, query_name, subs_mat, params):
    """ Perform pairwise local alignment using parsail-python """
    aln = parasail.sw_trace_striped_32(query, reference, params['gap_open'], params['gap_extend'], subs_mat)
    return process_alignment(aln, query, query_name, params)


def refine_locations(read, all_primers, locations, aln_params=DEFAULT_ALIGN_PARAMS, subs_mat=DEFAULT_SUBS_MAT):
    "Refine alignment edges based on local alignment"
    seq = read.Seq

    proc_locations = []
    for loc in locations:
        aln = pair_align(seq[loc.RefStart:loc.RefEnd], all_primers[loc.Query], loc.Query, subs_mat, aln_params)
        rscore = aln['norm_score']
        rloc = Hit(loc.Ref, loc.RefStart + aln["ref_start"], loc.RefStart + aln["ref_end"], loc.Query, aln["query_start"], aln["query_end"], rscore)
        proc_locations.append(rloc)
    return proc_locations

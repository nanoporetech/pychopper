# -*- coding: utf-8 -*-

import edlib

from pychopper.common_structures import Hit
from pychopper import utils
from pychopper.parasail_backend import refine_locations


def find_locations(reads, all_primers, max_ed, pool, min_batch):
    "Find alignment hits of all primers in all reads using the edlib/parasail backend"
    for batch in utils.batch(reads, min_batch):
        for res in pool.map(_find_locations_single, zip(batch, [(all_primers, max_ed)] * len(batch))):
            try:
                yield res
            except StopIteration:
                return


def _find_locations_single(params):
    "Find alignment hits of all primers in a single reads using the edlib/parasail backend"
    read = params[0]
    all_primers, max_ed = params[1]
    all_locations = []
    for primer_acc, primer_seq in all_primers.items():
        primer_max_ed = int(max_ed * len(primer_seq))
        result = edlib.align(primer_seq, read.Seq,
                             mode="HW", task="locations", k=primer_max_ed)
        ed = result["editDistance"]
        locations = result["locations"]
        if locations:
            # all_locations[primer_acc] = []
            for refstart, refend in locations:
                refend += 1
                # ('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd Score')
                hit = Hit(read.Name, refstart, refend, primer_acc, 0,
                          len(primer_seq),  ed / len(primer_seq))
                all_locations.append(hit)
    refined_locations = refine_locations(read, all_primers, all_locations)
    return refined_locations

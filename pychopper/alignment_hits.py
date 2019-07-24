# -*- coding: utf-8 -*-

# Place here the code for sorting hits and removing overlaps.

from pychopper import seq_utils
from pychopper import hmmer_backend


def process_hits(hits, max_score):
    res = []
    for hit in sorted(hits, key=lambda x: (x.RefStart, x.RefEnd)):
        if hit.Score > max_score:
            continue
        if len(res) == 0:
            res.append(hit)
        else:
            last = res[-1]
            if last.RefEnd > hit.RefStart and hit.Score < last.Score:
                # Swap with top element:
                res[-1] = hit
                pass
            else:
                res.append(hit)
    return tuple(res)


if __name__ == "__main__":
    for s in seq_utils.readfq(open("tmp/100r.fq", "r")):
        hits = hmmer_backend.find_locations(s, "tmp/primers.hmm", E=1.0)
        proc_hits = process_hits(hits, 0.1)
        for ph in proc_hits:
            print(ph.RefStart, ph.RefEnd)

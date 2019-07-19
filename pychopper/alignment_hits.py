# -*- coding: utf-8 -*-


def process_hits(hits, max_score):
    "Process alignment hits by removing overlaps"
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

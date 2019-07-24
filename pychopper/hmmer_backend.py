# -*- coding: utf-8 -*-

import sys
import re
import subprocess as sp
from pychopper.common_structures import Hit, Seq
from pychopper import seq_utils


def _parse_hmmscan_tab(lines):
    res = []
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        tmp = re.split('\s+', line)
        hit = Hit(tmp[2], int(tmp[6]), int(tmp[7]) + 1, tmp[0], int(tmp[4]), int(tmp[5]) + 1, float(tmp[12]))
        res.append(hit)
    return res


def find_locations(read, phmm_file, E=0.1, threads=1):
    """
    """
    cmd = "nhmmscan --notextw --max -E {} --cpu {} --watson -o /dev/null --tblout /dev/stdout {} -"
    cmd = cmd.format(E, threads, phmm_file)
    proc = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE)
    dout, derr = proc.communicate(">{}\n{}\n".format(read.Name, read.Seq).encode())

    if derr is not None:
        print("Failed to run hmmscan with model {} and read {}:".format(phmm_file, read.Name), file=sys.stderr)
        sys.exit(1)

    hits = _parse_hmmscan_tab(dout.decode().split("\n"))
    return hits


if __name__ == "__main__":
    for s in seq_utils.readfq(open("tmp/100r.fq", "r")):
        find_locations(s, "tmp/primers.hmm")


from collections import OrderedDict
import subprocess as sp
from itertools import islice, chain
import numpy as np
import sys


def parse_config_string(s):
    res = OrderedDict()
    s = s.strip()
    for token in s.split("|"):
        token = token.strip().split(":")
        if len(token) != 2:
            raise Exception("Invalid token: " + token)
        strand = token[0].strip()
        if strand not in ('+', '-'):
            raise Exception("Invalid strand: " + strand)
        token = token[1].strip().split(",")
        res[(token[0].strip(), token[1].strip())] = strand
    return res


def batch(iterable, size):
    sourceiter = iter(iterable)
    while True:
        batchiter = islice(sourceiter, size)
        try:
            yield list(chain([next(batchiter)], batchiter))
        except StopIteration:
            return


def hit2bed(hit, read):
    # Hit = namedtuple('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd
    # Score')
    max_q, q = 100, None
    if hit.Score == 0:
        q = max_q
    elif hit.Score > 1:
        q = 0
    else:
        q = -10 * np.log10(hit.Score)
    name = hit.Query
    strand = "+"
    if hit.Query[0] == "-":
        name = hit.Query[1:]
        strand = "-"
    bed_line = "%s\t%d\t%d\t%s\t%d\t%s" % (read.Name, hit.RefStart, hit.RefEnd, name, q, strand)
    return bed_line


def count_fastq_records(fname, size=128000000, opener=open):
    fh = opener(fname, "r")
    count = 0
    while True:
        b = fh.read(size)
        if not b:
            break
        count += b.count("\n+\n")
    fh.close()
    return count


def check_command(cmd):
    if sp.call(cmd, shell=True) != 0:
        sys.stderr.write("Required command {} not found in the path!\n".format(cmd.split(" ")[0]))


def check_min_hmmer_version(major, minor):
    ver_line = sp.check_output("nhmmscan -h | grep '^# HMMER'", shell=True).decode()
    tmp = ver_line.split("HMMER")[1].strip().split(" ")[0].split(".")
    f_major, f_minor = int(tmp[0]), int(tmp[1])
    if f_major < major:
        raise Exception("Insufficient major HMMER version: " + str(f_major))
    if f_minor < minor:
        raise Exception("Insufficient minor HMMER version: " + str(f_minor))

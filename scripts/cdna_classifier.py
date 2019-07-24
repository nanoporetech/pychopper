#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import tqdm
from pychopper import seq_utils as seu
from pychopper import utils
from pychopper import chopper
from pychopper import report
import pandas as pd
from collections import OrderedDict, defaultdict

"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(
    description='Tool to identify full length cDNA reads. Primers have to be specified as they are on the forward strand.')
parser.add_argument(
    '-b', metavar='primers', type=str, default=None, help="Primers fasta.", required=False)
parser.add_argument(
    '-g', metavar='phmm_file', type=str, default=None, help="File with custom profile HMMs (None).", required=False)
parser.add_argument(
    '-c', metavar='config_string', type=str, default="+:SSP,-VNP|-:VNP,-SSP", help="Specify primer configurations for each direction (None).")
parser.add_argument(
    '-q', metavar='cutoff', type=float, default=98, help="Cutoff parameter (98).")
parser.add_argument(
    '-z', metavar='identity_cutoff', type=float, default=0.8, help="Lower identity threshold for finding primer/barcode (0.8).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, default=None, help="Report PDF.")
parser.add_argument(
    '-u', metavar='unclass_output', type=str, default=None, help="Write unclassified reads to this file.")
parser.add_argument(
    '-w', metavar='rescue_output', type=str, default=None, help="Write rescued reads to this file.")
parser.add_argument(
    '-S', metavar='stats_output', type=str, default=None, help="Write statistics to this file.")
parser.add_argument(
    '-A', metavar='scores_output', type=str, default=None, help="Write alignment scores to this file.")
parser.add_argument(
    '-m', metavar='method', type=str, default="edlib", help="Detection method (edlib or phmm).")
parser.add_argument(
    '-t', metavar='threads', type=int, default=4, help="Number of threads to use (4).")
parser.add_argument('input_fastx', metavar='input_fastx', type=str, help="Input file.")
parser.add_argument('output_fastx', metavar='output_fastx', type=str, help="Output file.")


def _new_stats():
    st = OrderedDict()
    st["Classification"] = OrderedDict([('Classified', 0), ('Rescue', 0), ('Unclassified', 0)])
    st["Strand"] = OrderedDict([('+', 0), ('-', 0)])
    st["RescueStrand"] = OrderedDict([('+', 0), ('-', 0)])
    st["RescueSegmentNr"] = defaultdict(int)
    st["RescueHitNr"] = defaultdict(int)
    st["UnclassHitNr"] = defaultdict(int)
    st["PercentUsable"] = defaultdict(int)
    return st


def _update_stats(st, segments, hits, usable_len, read):
    if len(segments) == 0:
        st["Classification"]["Unclassified"] += 1
        st["UnclassHitNr"][len(hits)] += 1
    elif len(segments) == 1:
        st["Classification"]["Classified"] += 1
        st["Strand"][segments[0].Strand] += 1
        st["PercentUsable"][int(segments[0].Len / len(read.Seq) * 100)] += 1
    else:
        for rs in segments:
            st["Classification"]["Rescue"] += 1
            st["RescueStrand"][rs.Strand] += 1
            st["RescueHitNr"][len(hits)] += 1
        st["PercentUsable"][int(sum([s.Len for s in segments]) / len(read.Seq) * 100)] += 1


if __name__ == '__main__':
    args = parser.parse_args()
    config = utils.parse_config_string(args.c)

    in_fh = sys.stdin
    if args.input_fastx != '-':
        in_fh = open(args.input_fastx, "r")

    out_fh = sys.stdout
    if args.output_fastx != '-':
        out_fh = open(args.output_fastx, "w")

    st = _new_stats()
    input_size = None
    if args.input_fastx != "-":
        input_size = os.stat(args.input_fastx).st_size
    pbar = tqdm.tqdm(total=input_size)

    if args.m == "edlib":
        all_primers, primer_length = seu.get_primers(args.b)
        max_ed = int(round( (1.0 - args.q) * primer_length))

    for read in seu.readfq(in_fh):
        if args.m == "phmm":
            segments, hits, usable_len = chopper.chopper_phmm(read, args.g, config, args.q, args.t)
        elif args.m == "edlib":
            segments, hits, usable_len = chopper.chopper_edlib(read, all_primers, config, max_ed, args.q)
        else:
            raise
        _update_stats(st, segments, hits, usable_len, read)
        for trim_read in chopper.segments_to_reads(read, segments):
            seu.writefq(trim_read, out_fh)
        pbar.update(seu.record_size(read, 'fastq'))

    print(st)
    in_fh.close()
    out_fh.flush()
    out_fh.close()
    pbar.close()

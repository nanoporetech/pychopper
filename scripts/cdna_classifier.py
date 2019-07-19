#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO
import os
import tqdm
from pychopper import seq_utils as seu
from pychopper import chopper
from pychopper import report
import pandas as pd
from collections import OrderedDict

"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(
    description='Tool to identify full length cDNA reads. Primers have to be specified as they are on the forward strand.')
parser.add_argument(
    '-b', metavar='primers', type=str, default=None, help="Primers fasta.", required=True)
parser.add_argument(
    '-g', metavar='phmm_dir', type=str, default=None, help="Directory with custom profile HMMs (None).", required=True)
parser.add_argument(
    '-i', metavar='input_format', type=str, default='fastq', help="Input/output format (fastq).")
parser.add_argument('-g', metavar='aln_params', type=str,
                    help="Alignment parameters (match, mismatch,gap_open,gap_extend).", default="1,-1,2,1")
parser.add_argument(
    '-q', metavar='cutoff', type=float, default=98, help="Cutoff parameter (98).")
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
parser.add_argument('input_fastx', metavar='input_fastx', type=str, help="Input file.")
parser.add_argument('output_fastx', metavar='output_fastx', type=str, help="Output file.")


# FIXME: the utility function below should be moved to a different file with documentation.
# FIXME: we need a way to specify the valid configuration for a read in each direction!

def _revcomp_seq(read):
    """ Reverse complement sequence and fix Id. """
    rev_read = read.reverse_complement()
    rev_read.id, rev_read.description = read.id, read.description
    return rev_read


def _filter_and_annotate(read, match):
    """ Filter sequences by match and annotate with direction. """
    if match is not None:
        direction = '+' if match == 'fwd_match' else '-'
        read.description = read.description + " strand=" + direction
        if match == 'rev_match':
            read = _revcomp_seq(read)
        return read, True
    else:
        return read, False


def _get_runid(desc):
    """ Parse out runid from sequence description. """
    tmp = [t for t in desc.split(" ") if t.startswith("runid")]
    if len(tmp) != 1:
        return "NA"
    return tmp[0].rsplit("=", 1)[1]


def _parse_aln_params(pstr):
    """ Parse alignment parameters. """
    res = {}
    tmp = [int(x) for x in pstr.split(',')]
    res['match'] = tmp[0]
    res['mismatch'] = tmp[1]
    res['gap_open'] = tmp[2]
    res['gap_extend'] = tmp[3]
    return res


def _record_size(read, in_format):
    """ Calculate record size. """
    dl = len(read.description)
    sl = len(read.seq)
    if in_format == 'fastq':
        bl = dl + 2 * sl + 6
    elif in_format == 'fasta':
        bl = dl + sl + 3
    else:
        raise Exception("Unkonwn format!")
    return bl


if __name__ == '__main__':
    args = parser.parse_args()

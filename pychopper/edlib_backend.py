# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import os
import sys
import argparse

import edlib

# import seq_utils
from pychopper.common_structures import Hit, Seq
from pychopper import seq_utils as seu


def find_locations(read, all_primers, max_ed = 5):
    """
        Input: A reads sequence represented as a string and a dictionary with primers. The dictionary has the primer accession as key and sequence as value.
        returns: Dictionary with primer accession as key and a list of tuples of the following form (start, stop, edit distance) as value.
    """

    all_locations = []
    for primer_acc, primer_seq in all_primers.items():
        result = edlib.align(primer_seq, read.Seq, mode="HW", task="locations", k=max_ed)
        ed = result["editDistance"]
        locations =  result["locations"]
        if locations:
            # all_locations[primer_acc] = []
            for refstart, refend in locations:
                # all_locations[primer_acc].append( (start, stop, ed))
                # ('Hit', 'Ref RefStart RefEnd Query QueryStart QueryEnd Score')
                hit = Hit(read.Name, refstart, refend, primer_acc, 0, len(primer_seq),  ed/(refend - refstart - 1))
                all_locations.append( hit )

    return all_locations



'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
    Much faster parsing of large files compared to Biopyhton.
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def get_primers_rev_comp(primers):
    all_primers = {}
    for acc, (seq, _) in readfq(open(args.primers, 'r')):
        all_primers[acc] = seq
        all_primers[ '-' + acc ] = seu.reverse_complement(seq)
    return all_primers, len(seq)


def main(args):
    all_primers, primer_length = get_primers_rev_comp(args.primers)
    # primer_length = len(all_primers.values()[0]) 
    # k = args.k
    k = int(round( (1.0 - args.q) * primer_length))
    print("edit distance set to:", k)
    #sys.exit()
    tot_no_full_length = 0
    for i, read in enumerate( seu.readfq(open(args.fastq, 'r'))):
        all_locations = find_locations(read, all_primers, k)
        if all_locations:
            has_1, has_2 = False, False
            for h in all_locations:
                if 'cDNA|1' in h.Query:
                    has_1 = True
                if 'cDNA|2' in h.Query:
                    has_2 = True
            if not (has_1 and has_2):
                tot_no_full_length += 1
                # print("read {0} did not have both barcodes".format(i), all_locations)
        if i % 100000 == 0:
            print("tot iterated: {0}, tot not having both primers:{1}".format(i, tot_no_full_length))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="De novo clustering of long-read transcriptome reads", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fastq', type=str,  default=False, help='Path to input fastq folder with reads in clusters')
    parser.add_argument('--primers', type=str,  default=False, help='Path to input fastq folder with primers')
    parser.add_argument('--k', type=int, default=5, help='Max edit distnace')
    parser.add_argument('--q', type=float, default=0.8, help='Minimum identity')
    args = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    
    main(args)


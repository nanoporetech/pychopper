# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import os
import sys
import argparse

import edlib

import seq_utils



def find_locations(read, all_primers, max_ed = 5):
    """
        Input: A reads sequence represented as a string and a dictionary with primers. The dictionary has the primer accession as key and sequence as value.
        returns: Dictionary with primer accession as key and a list of tuples of the following form (start, stop, edit distance) as value.
    """

    all_locations = {}
    for acc, primer_seq in all_primers.items():
        result = edlib.align(primer_seq, read, mode="HW", task="locations", k=max_ed)
        ed = result["editDistance"]
        locations =  result["locations"]
        if locations:
            all_locations[acc] = []
            for start, stop in locations:
                all_locations[acc].append( (start, stop, ed))

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
        all_primers[acc + '_rc'] = seq_utils.reverse_complement(seq)
    return all_primers


def main(args):
    all_primers = get_primers_rev_comp(args.primers)
    k = args.k
    for i, (acc, (seq, qual)) in enumerate(readfq(open(args.fastq, 'r'))):
        all_locations = find_locations(seq, all_primers, k)
        if all_locations:
            print("read {0} had barcode".format(i), all_locations)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="De novo clustering of long-read transcriptome reads", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fastq', type=str,  default=False, help='Path to input fastq folder with reads in clusters')
    parser.add_argument('--primers', type=str,  default=False, help='Path to input fastq folder with primers')
    parser.add_argument('--k', type=int, default=5, help='Max edit distnace')
    args = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    
    main(args)


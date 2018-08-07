# -*- coding: utf-8 -*-

from six.moves import reduce
import sys
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

""" Utilities manipulating biological sequences and formats. Extensions to biopython functionality.
"""

# Reverse complements of bases, taken from dragonet:
comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    '-': '-'
}

# Shortcut to the list of DNA bases:
bases = sorted(list(IUPACUnambiguousDNA().letters))
ambiguous_bases = sorted(list(IUPACAmbiguousDNA().letters))


def base_complement(k):
    """ Return complement of base.

    Performs the subsitutions: A<=>T, C<=>G, X=>X for both upper and lower
    case. The return value is identical to the argument for all other values.

    :param k: A base.
    :returns: Complement of base.
    :rtype: str

    """
    try:
        return comp[k]
    except KeyError:
        sys.stderr.write(
            "WARNING: No reverse complement for {} found, returning argument.".format(k))
        return k


def reverse_complement(seq):
    """ Return reverse complement of a string (base) sequence.

    :param seq: Input sequence.
    :returns: Reverse complement of input sequence.
    :rtype: str

    """
    if len(seq) == 0:
        return seq
    return reduce(lambda x, y: x + y, map(base_complement, seq[::-1]))


def mock_qualities(record, mock_qual):
    """Add mock quality values to SeqRecord object.

    :param record: A SeqRecord object.
    :param mock_qual: Mock quality value used for each base.
    :returns: The record augmented with mock quality values.
    :rtype: object

    """
    rec_copy = record[:]
    rec_copy.letter_annotations["phred_quality"] = [mock_qual] * len(rec_copy)
    return rec_copy


def new_dna_record(sequence, name, qualities=None):
    """Create a new SeqRecord object using IUPACUnambiguousDNA and the specified sequence.

    :param sequence: The sequence.
    :param name: Record identifier.
    :param qualities: List of base qualities.
    :returns: The SeqRecord object.
    :rtype: SeqRecord

    """
    seq_record = SeqRecord(
        Seq(sequence, IUPACUnambiguousDNA), id=name, description="", name="")
    if qualities is not None:
        seq_record.letter_annotations["phred_quality"] = qualities
    return seq_record


def write_seq_records(records_iterator, output_object, format='fasta'):
    """Write out SeqRecord objects to a file from an iterator in the specified format.

    :param records_iterator: An iterator of SeqRecord objects.
    :param output_object: Open file object or file name.
    :param format: Output format (fasta by default).
    :returns: None
    :rtype: object

    """
    if type(output_object) != str:
        SeqIO.write(records_iterator, output_object, format)
    else:
        with open(output_object, 'w') as output_handle:
            SeqIO.write(records_iterator, output_handle, format)


def read_seq_records(input_object, format='fasta'):
    """Read SeqRecord objects from a file in the specified format.

    :param input_object: A file object or a file name.
    :param format: Input format (fasta by default).
    :returns: A dictionary with the parsed SeqRecord objects.
    :rtype: generator

    """
    handle = input_object
    if type(handle) == str:
        handle = open(handle, "rU")
    return SeqIO.parse(handle, format)


def count_records(input_object, format='fasta'):
    """Count SeqRecord objects from a file in the specified format.

    :param input_object: A file object or a file name.
    :param format: Input format (fasta by default).
    :returns: Number of records in input file.
    :rtype: int

    """
    handle = input_object
    if type(handle) == str:
        handle = open(handle, "rU")
    counter = 0
    for _ in SeqIO.parse(handle, format):
        counter += 1
    return counter


def read_seq_records_dict(input_object, format='fasta'):
    """Read SeqRecord objects to a dictionary from a file in the specified format.

    :param input_object: A file object or a file name.
    :param format: Input format (fasta by default).
    :returns: An iterator of SeqRecord objects.
    :rtype: dict

    """
    handle = input_object
    if type(handle) == str:
        handle = open(handle, "rU")
    return SeqIO.to_dict(SeqIO.parse(handle, format))

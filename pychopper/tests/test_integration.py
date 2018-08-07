import unittest
import os
from os import path
import subprocess


class TestIntegration(unittest.TestCase):

    def testIntegration(self):
        """ Integration test. """
        base = path.dirname(__file__)
        pkg_base = path.dirname(path.dirname(base))
        test_base = path.join(base, 'data')

        cdna_classifier = path.join(pkg_base, 'scripts', 'cdna_classifier.py')
        barcodes = path.join(test_base, 'barcodes.fas')
        input_fasta = path.join(test_base, 'ref.fas')
        output_fasta = path.join(test_base, 'test_output.fas')
        expected_output = path.join(test_base, 'expected_output.fas')

        subprocess.call([cdna_classifier, '-i', 'fasta','-s','95.0', '-b', barcodes, input_fasta, output_fasta])
        retval = subprocess.call(['cmp', output_fasta, expected_output])
        self.assertEqual(retval, 0)
        os.remove(output_fasta)

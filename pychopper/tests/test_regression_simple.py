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

        cdna_classifier = path.join(pkg_base, 'scripts', 'pychopper')
        barcodes = path.join(test_base, 'barcodes.fas')
        input_fasta = path.join(test_base, 'ref.fq')
        output_fasta = path.join(test_base, 'test_output.fq')
        expected_output = path.join(test_base, 'expected_output.fas')

        subprocess.call("{} {} {} {} {}".format(cdna_classifier, "-Y 0 -B 3 -q 0.5 -m edlib -b", barcodes, input_fasta, output_fasta), shell=True)
        retval = subprocess.call(['cmp', output_fasta, expected_output])
        self.assertEqual(retval, 0)
        os.remove(output_fasta)

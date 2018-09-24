Pychopper: A tool to identify full length cDNA reads
====================================================

Installation
------------

Install via pip:

```
pip install git+https://github.com/nanoporetech/pychopper.git
```

Or clone the repository:

```
git clone https://github.com/nanoporetech/pychopper.git
```

And install the package:

```
python setup.py install
```

Install the package in developer mode:

```
python setup.py develop
```

Run the tests:

```
make test
```

Build the documentation:

```
make docs
```

Issue `make help` to get a list of `make` targets.

Usage
-----

```
usage: cdna_classifier.py [-h] -b barcodes [-i input_format] [-g aln_params]
                          [-t target_length] [-s score_percentile]
                          [-n sample_size] [-r report_pdf] [-u unclass_output]
                          input_fastx output_fastx

Tool to identify full length cDNA reads. Primers have to specified as they are
on the forward strand.

positional arguments:
  input_fastx          Input file.
  output_fastx         Output file.

optional arguments:
  -h, --help           show this help message and exit
  -b barcodes          Primers fasta.
  -i input_format      Input/output format (fastq).
  -g aln_params        Alignment parameters (match,
                       mismatch,gap_open,gap_extend).
  -t target_length     Number of bases to scan at each end (200).
  -s score_percentile  Score cutoff percentile (98).
  -n sample_size       Number of samples when calculating score cutoff
                       (100000).
  -r report_pdf        Report PDF.
  -u unclass_output    Write unclassified reads to this file.
```

Example usage:

```bash
cdna_classifier.py -b cdna_barcodes.fas -r report.pdf -u unclassified.fq input.fq full_length_output.fq
```

The primers have to specified as they are on the forward strand (see `data/cdna_barcodes.fas` for an example).
The score cutoffs for each primer are calculated by aligning them against random sequences and taking the `-s` percentile of the score distribution (98 by default).

Contributing
------------

- Please fork the repository and create a merge request to contribute.
- Use [bumpversion](http://bit.ly/2cSUryt) to manage package versioning.
- The code should be [PEP8](https://www.python.org/dev/peps/pep-0008) compliant, which can be tested by `make lint`.

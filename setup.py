#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from glob import glob
import re
import os
import pkg_resources

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'edlib',
    'parasail',
    'matplotlib',
    'seaborn',
    'tqdm',
    'six',
    'pandas',
    'pytest',
    'sphinx',
    'sphinx_rtd_theme',
]

test_requirements = [
    # TODO: put package test requirements here
]

__pkg_name__ = 'pychopper'
__author__ = 'cwright'
__description__ = 'Identify, orient, and trim full-length ONT cDNA reads.'

__path__ = os.path.dirname(__file__)
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_name__))
# Get the version number from __init__.py
verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)

if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))

dir_path = os.path.dirname(__file__)
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    install_requires = [
        str(requirement) for requirement in
        pkg_resources.parse_requirements(fh)]

data_files = []
extra_requires = {}
extensions = []

setup(
    name='pychopper',
    version= __version__,
    description="A tool to identify full length cDNA reads.",
    long_description=readme,
    author="ONT Applications Group",
    author_email='Apps@nanoporetech.com',
    url='',
    packages=[
        'pychopper',
        'pychopper.phmm_data',
        'pychopper.primer_data'
    ],
    package_dir={'pychopper':
                 'pychopper'},
    package_data={'pychopper': ['primer_data/*.fas', 'phmm_data/*.*']},
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='pychopper',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
    ],
    tests_require=test_requirements,
    scripts=[x for x in glob('scripts/*.py') if x != 'scripts/__init__.py']
)

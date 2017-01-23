#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from collections import OrderedDict
import itertools
import sys
sys.path.append('..')

from uncle_PSL import psl2sam

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Script to convert PSL files (BLAT output) to SAM format.')
parser.add_argument(
    '-q', metavar='fastq', type=str, help="Input fastq.", required=False, default=None)
parser.add_argument(
    'infile', metavar='input_file', type=str, help="Input file.")

# Reference on the PSL format: http://www.ensembl.org/info/website/upload/psl.html
# Reference on the SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf


if __name__ == '__main__':
    args = parser.parse_args()

    psl2sam.psl2sam(args.infile)

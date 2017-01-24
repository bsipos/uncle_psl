#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# (c) 2016 Oxford Nanopore Technologies Ltd.

import argparse
from collections import OrderedDict
import itertools
import sys

from Bio import SeqIO

from uncle_PSL import psl2sam

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Script to convert PSL files (BLAT output) to SAM format.')
parser.add_argument(
    '-f', metavar='reads_fasta', type=str, help="Reads in fasta format.", required=False, default=None)
parser.add_argument(
    '-N', metavar='n_limit', type=int, help="Use N CIGAR operation for deletions larger than this parameter (None).", required=False, default=None)
parser.add_argument(
    '-H', action="store_false", help="Use hard clipping instead of soft clipping.", default=True)
parser.add_argument('infile', nargs='?', help='Input PSL (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('outfile', nargs='?', help='Output SAM (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)

# Reference on the PSL format: http://www.ensembl.org/info/website/upload/psl.html
# Reference on the SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf


if __name__ == '__main__':
    args = parser.parse_args()
    reads = SeqIO.to_dict(SeqIO.parse(args.f, 'fasta')) if args.f is not None else None
    psl2sam.psl2sam(args.infile, args.outfile, reads, args.H, args.N)

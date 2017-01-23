# -*- coding: utf-8 -*-

from collections import OrderedDict
import itertools

def _prepare_psl_dict():
    fields_text = """matches
misMatches
repMatches
nCount
qNumInsert
qBaseInsert
tNumInsert
tBaseInsert
strand
qName
qSize
qStart
qEnd
tName
tSize
tStart
tEnd
blockCount
blockSizes
qStarts
tStarts"""
    fields = fields_text.split("\n")
    fields_dict = OrderedDict()
    for field in fields:
        fields_dict[field] = None
    return fields_dict


def _iter_fields(fname, nr_fields=21):
    fh = open(fname, 'r')
    for line in fh:
        fields = line.split()
        if len(fields) != nr_fields:
            continue
        yield fields


def _generate_cigar(qStart, blockSizes, qStarts, tStarts, blockCount, qSize, qEnd):
    # Construct the CIGAR string, here we go:
    cigar = []

    # 5' hard clipping:
    if qStart != 0:
        cigar.append("{}H".format(qStart))
    # Process alignment segments:
    bs = blockSizes[0]
    qs = qStarts[0]
    ts = tStarts[0]
    for i in xrange(1, blockCount):  # process block by block, start from 2nd block
        # Gap between query blocks:
        insertion = qStarts[i] - qStarts[i - 1] - bs
        # Gap between target blocks:
        deletion = tStarts[i] - tStarts[i - 1] - bs
        # Add block to cigar as match:
        cigar.append("{}M".format(bs))
        # Add inserion:
        if insertion != 0:
            cigar.append("{}I".format(insertion))
        # Add deletion:
        if deletion != 0:
            cigar.append("{}D".format(deletion))
        # NM?
        # Advance to next block:
        bs, qs, ts = blockSizes[i], qStarts[i], tStarts[i]
    # Add last block:
    cigar.append("{}M".format(bs))
    # Deal with 3' clipping:
    three_clip = qSize - qEnd
    if three_clip != 0:
        cigar.append("{}H".format(three_clip))
    cigar = ''.join(cigar)
    return cigar


def psl_rec2sam_rec(psl):
    # Figure out strand:
    if len(psl['strand']) != 2:
        raise Exception('Invalid strand field in record: {}'.format(ps['qName']))
    strand = '+' if all(x == list(psl['strand'])[0] for x in list(psl['strand'])) else '-'

    # Type conversion of coordinates:
    qStart, qEnd = int(psl['qStart']), int(psl['qEnd'])
    tStart, tEnd = int(psl['tStart']), int(psl['tEnd'])
    qSize, tSize = int(psl['qSize']), int(psl['tSize'])

    # Transform start and end coordinates if strand is '-':
    if strand == '-':
        qStart = qSize - qEnd
        qEnd = qSize - int(psl['qStart'])

    # Extract segement information:
    blockCount = int(psl['blockCount'])
    blockSizes = [int(bs) for bs in psl['blockSizes'].split(',') if len(bs) > 0]
    qStarts = [int(qs) for qs in psl['qStarts'].split(',') if len(qs) > 0]
    tStarts = [int(ts) for ts in psl['tStarts'].split(',') if len(ts) > 0]

    # Reverse and transform segment information if strand is '-':
    if psl['strand'][1] == '-':
        blockSizes = blockSizes[::-1]
        qStarts = qStarts[::-1]
        tStarts = tStarts[::-1]
        for i in xrange(blockCount):
            qStarts[i] = qSize - blockSizes[i] - qStarts[i]
            tStarts[i] = tSize - blockSizes[i] - tStarts[i]

    # Generate CIGAR:
    cigar = _generate_cigar(qStart, blockSizes, qStarts, tStarts, blockCount, qSize, qEnd)
    print cigar

def psl2sam(psl_file):
    for fields in _iter_fields(psl_file):
        psl_fields = _prepare_psl_dict()
        for pos, key in enumerate(psl_fields.keys()):
            psl_fields[key] = fields[pos]
        psl_rec2sam_rec(psl_fields)

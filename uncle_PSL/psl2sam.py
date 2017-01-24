# -*- coding: utf-8 -*-


# Reference on the PSL format: http://www.ensembl.org/info/website/upload/psl.html
# Reference on the SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf

from collections import OrderedDict
import itertools

from uncle_PSL.sam_writer import SamWriter
from uncle_PSL.seq_util import reverse_complement


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


def _iter_fields(handle, nr_fields=21):
    for line in handle:
        fields = line.split()
        if len(fields) != nr_fields:
            continue
        yield fields


def _generate_cigar(qStart, blockSizes, qStarts, tStarts, blockCount, qSize, qEnd, strand, soft_clip=True, n_limit=None):
    # Construct the CIGAR string, here we go:
    cigar = []
    indels = 0
    clip_op = 'S' if soft_clip else 'H'

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
            del_op = 'D'
            # Use N operation if deleltion is larger than limit:
            if (n_limit is not None) and (deletion >= n_limit):
                del_op = 'N'
            cigar.append("{}{}".format(deletion, del_op))
        indels += deletion
        indels += insertion
        # Advance to next block:
        bs, qs, ts = blockSizes[i], qStarts[i], tStarts[i]
    # Add last block:
    cigar.append("{}M".format(bs))
    # 5' hard clipping:
    if qStart != 0:
        cigar.insert(0,"{}{}".format(qStart, clip_op))
    # Deal with 3' clipping:
    three_clip = qSize - qEnd
    if three_clip != 0:
        cigar.append("{}{}".format(three_clip, clip_op))
    # CIGAR complete:
    return cigar, indels

def _extract_segment_info(psl):
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

    return blockCount, blockSizes, qStarts, tStarts


def psl_rec2sam_rec(psl, sam_writer, reads, soft_clip, n_limit):
    # Figure out strand:
    if len(psl['strand']) == 1:
        strand = psl['strand']  # Not sure if this is sane!
        psl['strand'] = psl['strand'] + '+' # Assume +
    if len(psl['strand']) == 2:
        strand = '+' if all(x == list(psl['strand'])[0] for x in list(psl['strand'])) else '-'
    else:
        raise Exception('Invalid strand field in record: {}'.format(psl['qName']))

    # Type conversion of coordinates:
    qStart, qEnd = int(psl['qStart']), int(psl['qEnd'])
    tStart, tEnd = int(psl['tStart']), int(psl['tEnd'])
    qSize, tSize = int(psl['qSize']), int(psl['tSize'])

    # Transform start and end coordinates if strand is '-':
    if strand == '-':
        qStart = qSize - qEnd
        qEnd = qSize - int(psl['qStart'])

    # Extract segement information:
    blockCount, blockSizes, qStarts, tStarts = _extract_segment_info(psl)

    # Generate CIGAR:
    cigar, indels = _generate_cigar(
        qStart, blockSizes, qStarts, tStarts, blockCount, qSize, qEnd, strand, soft_clip, n_limit)
    cigar_string = ''.join(cigar)
    NM = indels + int(psl['misMatches']) + int(psl['nCount'])

    # Construct SAM record:
    flag = 0 if strand == '+' else 16  # Strand flag
    # Construct sequence:
    seq = '*'
    if reads is not None and psl['qName'] in reads:
        seq = str(reads[psl['qName']].seq)
        if strand == '-':
            seq = reverse_complement(seq)
    # Deal with hard clipping (code could be cleaner):
    # Clip 5':
    first_op = cigar[0]
    if first_op[-1] == 'H':
        seq = seq[int(first_op[:-1]):]
    last_op = cigar[-1]
    # Clip 3':
    if last_op[-1] == 'H':
        seq = seq[:len(seq) - int(last_op[:-1])]

    sam = sam_writer.new_sam_record(qname=psl['qName'], flag=flag, rname=psl['tName'], pos=int(psl['tStart']) + 1,
                                    mapq=0, cigar=cigar_string, rnext='*', pnext=0, tlen=0, seq=seq, qual='*', tags='NM:i:{}'.format(NM))
    return sam


def psl2sam(psl_handle, out_handle, reads, soft_clip=True, n_limit=None):
    # Create SamWriter object:
    sam_writer = SamWriter(out_handle)
    # Iterate PSL records:
    for fields in _iter_fields(psl_handle):
        psl_fields = _prepare_psl_dict()
        # Fill PSL structure:
        for pos, key in enumerate(psl_fields.keys()):
            psl_fields[key] = fields[pos]
        # Convert PSL -> SAM:
        sam_rec = psl_rec2sam_rec(psl_fields, sam_writer, reads, soft_clip, n_limit)
        sam_writer.write(sam_rec)

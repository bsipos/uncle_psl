import unittest
from os import path
import tempfile

from uncle_PSL import psl2sam


class ExamplePsl2sam(unittest.TestCase):

    def _parse_sam(self, fname):
        """ Simple parsing of SAM files. """
        fh = open(fname, 'r')
        lines = [l.split() for l in fh if not l.startswith('@')]
        res = []
        for l in lines:
            res.append({'ref': l[0], 'flag': l[1], 'query': l[2], 'pos': l[3], 'cigar': l[5]})
        return res

    def test_psl2sam(self):
        """ Test PSL -> SAM conversion against BWA results. """
        top = path.dirname(__file__)
        bwa_sam = path.join(top, "data/bwa.sam")
        psl = path.join(top, "data/blat_top.psl")
        res_sam = tempfile.NamedTemporaryFile(prefix='test_psl2sam', delete=False)
        psl_fh = open(psl, 'r')
        psl2sam.psl2sam(psl_fh, res_sam, None, soft_clip=True, n_limit=None)
        res_sam.flush()
        bwa_records = self._parse_sam(bwa_sam)
        psl_records = self._parse_sam(res_sam.name)
        res_sam.close()
        self.assertEqual(bwa_records, psl_records)

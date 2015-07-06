"""
Tests the parsers in the parsers module
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

import ga4gh.parsers.cufflinks_parser as cufflinks_parser
import ga4gh.parsers.rsem_parser as rsem_parser


class TestCufflinksParser(unittest.TestCase):

    def testParser(self):
        # check that we can parse an example file without errors
        path = 'tests/data/cufflinks.txt'
        with file(path) as dataFile:
            parser = cufflinks_parser.CufflinksParser(dataFile)
            rows = parser.parse()
        self.assertEqual(len(rows), 49)

        # test parsed values of an arbitrary row
        row = rows[2]
        self.assertEqual(row.tracking_id, 'ENSMUSG00000051951')
        self.assertEqual(row.class_code, '-')
        self.assertEqual(row.nearest_ref_id, '-')
        self.assertEqual(row.gene_id, 'ENSMUSG00000051951')
        self.assertEqual(row.gene_short_name, 'Xkr4')
        self.assertEqual(row.tss_id, '-')
        self.assertEqual(row.locus, 'chr1:3195981-3661579')
        self.assertEqual(row.length, '-')
        self.assertEqual(row.coverage, '-')
        self.assertAlmostEqual(row.FPKM, 0.391061)
        self.assertEqual(row.FPKM_conf_lo, 0.0546014)
        self.assertEqual(row.FPKM_conf_hi, 0.727521)
        self.assertEqual(row.FPKM_status, 'OK')


class TestRsemParser(unittest.TestCase):

    def testParser(self):
        # check that we can parse an example file without errors
        path = 'tests/data/abundance.trimmed.txt'
        with file(path) as dataFile:
            parser = rsem_parser.RsemParser(dataFile)
            rows = parser.parse()
        self.assertEqual(len(rows), 10)

        # test parsed values of an arbitrary row
        row = rows[0]
        self.assertEqual(row.target_id, "ENST00000390435")
        self.assertEqual(row.length, 341)
        self.assertEqual(row.eff_length, 162)
        self.assertEqual(row.est_counts, 1)
        self.assertAlmostEqual(row.tpm, 0.20407)

"""
A parser for the cufflinks data format
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.parsers.basic_parser as basic_parser


class CufflinksRow(basic_parser.AbstractRowObject):
    """
    A row in a cufflinks file
    """
    metainfo = (
        ('tracking_id', str),
        ('class_code', str),
        ('nearest_ref_id', str),
        ('gene_id', str),
        ('gene_short_name', str),
        ('tss_id', str),
        ('locus', str),
        ('length', str),
        ('coverage', str),
        ('FPKM', float),
        ('FPKM_conf_lo', float),
        ('FPKM_conf_hi', float),
        ('FPKM_status', str),
    )

    def __init__(self, values):
        super(CufflinksRow, self).__init__(values)


class CufflinksParser(basic_parser.BasicParser):
    """
    A parser for the cufflinks format
    """
    def processHeader(self):
        # skip the 1 row header
        self._stream.readline()

    def getRowClass(self):
        return CufflinksRow

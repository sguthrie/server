"""
A parser for the RSEM data format
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import ga4gh.parsers.basic_parser as basic_parser


class RsemRow(basic_parser.AbstractRowObject):
    """
    A row in a RSEM file
    """
    metainfo = (
        ('target_id', str),
        ('length', int),
        ('eff_length', int),
        ('est_counts', float),
        ('tpm', float),
    )

    def __init__(self, values):
        super(RsemRow, self).__init__(values)


class RsemParser(basic_parser.BasicParser):
    """
    A parser for the RSEM format
    """
    def processHeader(self):
        # skip the 1 row header
        self._stream.readline()

    def getRowClass(self):
        return RsemRow

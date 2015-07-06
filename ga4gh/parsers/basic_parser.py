"""
A basic parser for parsing bioinformatics data
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals


class ParseException(Exception):
    """
    Something went wrong during parsing
    """


class AbstractRowObject(object):
    """
    An abstract base class of a row object for a file that a subclass
    of BasicParser can parse.

    Subclasses must define a class-level 'metainfo' attribute that is a
    tuple of 2-tuples.  Each tuple is the name and class of an attribute
    corresponding to the row position equal to the metainfo list index.
    """
    def __init__(self, values):
        if len(values) != len(self.metainfo):
            msg = "Row has wrong number of values; has {}, should have {}"
            msgStr = msg.format(len(values), len(self.metainfo))
            raise ParseException(msgStr)
        for position, value in enumerate(values):
            attrMeta = self.metainfo[position]
            attrName = attrMeta[0]
            attrValueClass = attrMeta[1]
            attrValue = attrValueClass(value)
            setattr(self, attrName, attrValue)

    def toDict(self):
        """
        Returns a dictionary representation of the row
        """
        dic = {}
        for attrName, _ in self.metainfo:
            dic[attrName] = getattr(self, attrName)
        return dic


class BasicParser(object):
    """
    An abstract parser for parsing bioinformatics file formats.

    Subclasses must define processHeader and getRowClass methods.
    """
    def __init__(self, stream):
        self._stream = stream
        self._objects = []

    def processHeader(self):
        """
        Handles the processing of the file header
        """
        raise NotImplementedError()

    def parse(self):
        """
        Parses the file data.  Returns a list of row objects.
        """
        self.processHeader()
        rowClass = self.getRowClass()
        objects = []
        for line in self._stream:
            values = line.split()
            obj = rowClass(values)
            objects.append(obj)
        return objects

    def getRowClass(self):
        """
        Returns the class that should be used to create row objects.
        """
        raise NotImplementedError()

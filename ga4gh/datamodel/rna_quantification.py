"""
Module responsible for translating feature expression data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.protocol as protocol

"""
TODO: Would be nice to just use the csv module to read inputs and have headers in files for clarity
      and to eliminate the whole record[N] absurdity.
"""
class RNASeqResult(object):
    """
    Class representing a single RnaQuantification in the GA4GH data model.
    """
    def __init__(self, rnaQuantificationId, rnaQuantDataPath):
        self._rnaQuantificationId = rnaQuantificationId
        self._rnaQuantificationFile = os.path.join(rnaQuantDataPath, "rnaseq.table")
        self._characterizationFile = os.path.join(rnaQuantDataPath, "dist.table")
        self._readCountFile = os.path.join(rnaQuantDataPath, "counts.table")

    def convertCharacterization(self, record):
        readCharacterization = protocol.Characterization
        readCharacterization.analysisId = record[0]
        readCharacterization.complexity = float(record[1])
        readCharacterization.exonicFraction = float(record[2])
        readCharacterization.fractionMapped = float(record[3])
        readCharacterization.intergenicFraction = float(record[4])
        readCharacterization.intronicFraction = float(record[5])

        return readCharacterization

    def getCharacterization(self, rnaQuantificationId):
        """
        input is tab file with no header.  Columns are:
        analysisId, complexity, exonicFraction, fractionMapped, intergenicFraction, intronicFraction
        """
        characterizationData = open(self._characterizationFile, "r")
        quantCharacterization = characterizationData.readline()
        fields = quantCharacterization.split('/t')
        if rnaQuantificationID is None or fields[0] == rnaQuantificationId:
            yield self.convertCharacterization(fields)

    def convertReadCounts(self, record):
        readCount = protocol.ReadCounts
        readCount.analysisId = record[0]
        readCount.multiCount = int(record[1])
        readCount.multiSpliceCount = int(record[2])
        readCount.totalReadCount = int(record[3])
        readCount.uniqueCount = int(record[4])
        readCount.uniqueSpliceCount = int(record[5])

        return readCount

    def getReadCounts(self, rnaQuantificationId):
        """
        input is tab file with no header.  Columns are:
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount, uniqueSpliceCount
        """
        readCountData = open(self._readCountFile, "r")
        countData = readCountData.readline()
        fields = countData.split('/t')
        if rnaQuantificationID is None or fields[0] == rnaQuantificationId:
            yield self.convertReadCounts(fields)

    def convertRnaQuantification(self, record):
        rnaQuantification = protocol.RnaQuantification
        rnaQuantification.id = record[0]
        rnaQuantification.annotationIds = record[1].split(',')
        rnaQuantification.description = record[2]
        rnaQuantification.name = record[3]
        rnaQuantification.readGroupId = record[4]

        return rnaQuantification

    def getRnaQuantification(self, rnaQuantificationId):
        """
        input is tab file with no header.  Columns are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        rnaQuantificationData = open(self._rnaQuantificationFile, "r")
        quantData = rnaQuantificationData.readline()
        fields = quantData.strip().split('\t')
        if rnaQuantificationId is None or fields[0] == rnaQuantificationId:
            yield self.convertRnaQuantification(fields)


class SimulatedRNASeqResult(object):
    """
    An RNA Quantification that doesn't derive from a data store.
    Used mostly for testing.
    """
    def __init__(self, rnaQuantificationId, rnaQuantDataPath):
        self._rnaQuantificationId = rnaQuantificationId

    def generateCharacterization(self):
        """
            Currently just returns default values.
        """
        readCharacterization = protocol.Characterization

        return readCharacterization

    def getCharacterization(self, rnaQuantificationId):
        """
        input is tab file with no header.  Columns are:
        analysisId, complexity, exonicFraction, fractionMapped, intergenicFraction, intronicFraction
        """
        characterizationData = open(self._characterizationFile, "r")
        quantCharacterization = characterizationData.readline()
        fields = quantCharacterization.split('/t')
        if rnaQuantificationID is None or fields[0] == rnaQuantificationId:
            yield self.generateCharacterization()

    def generateReadCounts(self):
        """
            Currently just returns default values.
        """
        readCount = protocol.ReadCounts

        return readCount

    def getReadCounts(self, rnaQuantificationId):
        """
        input is tab file with no header.  Columns are:
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount, uniqueSpliceCount
        """
        readCountData = open(self._readCountFile, "r")
        countData = readCountData.readline()
        fields = countData.split('/t')
        if rnaQuantificationID is None or fields[0] == rnaQuantificationId:
            yield self.generateReadCounts()

    def generateRnaQuantification(self):
        """
            Currently just returns default values.
        """
        rnaQuantification = protocol.RnaQuantification

        return rnaQuantification

    def getRnaQuantification(self, rnaQuantificationId):
        """
        input is tab file with no header.  Columns are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        rnaQuantificationData = open(self._rnaQuantificationFile, "r")
        quantData = rnaQuantificationData.readline()
        fields = quantData.strip().split('\t')
        if rnaQuantificationId is None or fields[0] == rnaQuantificationId:
            yield self.generateRnaQuantification()

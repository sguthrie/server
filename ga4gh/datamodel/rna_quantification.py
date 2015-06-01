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
TODO: Would be nice to just use the csv module to read inputs and have headers
in files for clarity and to eliminate the whole record[N] absurdity.

Additionally, characterization and read counts doesn't have an access point.
"""


class RNASeqResult(object):
    """
    Class representing a single RnaQuantification in the GA4GH data model.
    """
    def __init__(self, rnaQuantificationId, rnaQuantDataPath):
        self._rnaQuantificationId = rnaQuantificationId
        self._rnaQuantificationFile = os.path.join(
            rnaQuantDataPath, "rnaseq.table")
        self._characterizationFile = os.path.join(
            rnaQuantDataPath, "dist.table")
        self._readCountFile = os.path.join(rnaQuantDataPath, "counts.table")
        self._expressionLevelFile = os.path.join(
            rnaQuantDataPath, "expression.table")

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
        analysisId, complexity, exonicFraction, fractionMapped,
        intergenicFraction, intronicFraction
        """
        characterizationData = open(self._characterizationFile, "r")
        quantCharacterization = characterizationData.readline()
        fields = quantCharacterization.split('/t')
        if rnaQuantificationId is None or fields[0] == rnaQuantificationId:
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
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount,
        uniqueSpliceCount
        """
        readCountData = open(self._readCountFile, "r")
        countData = readCountData.readline()
        fields = countData.split('/t')
        if rnaQuantificationId is None or fields[0] == rnaQuantificationId:
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

    def convertExpressionLevel(self, record):
        expressionLevel = protocol.ExpressionLevel()
        expressionLevel.id = record[0]
        expressionLevel.annotationId = record[1]
        expressionLevel.expression = record[2]
        expressionLevel.featureGroupId = record[3]
        expressionLevel.isNormalized = record[4]
        expressionLevel.rawReadCount = record[5]
        expressionLevel.score = record[6]
        expressionLevel.units = record[7]

        return expressionLevel

    def getExpressionLevel(self, expressionLevelId, featureGroupId):
        """
        input is tab file with no header.  Columns are:
        id, annotationId, expression, featureGroupId,
        isNormalized, rawReadCount, score, units

        expressionLevelId is not None: return only the specific expressionLevel
        object
        featureGroupId is not None: return all in that group
        """
        expressionLevelData = open(self._expressionLevelFile, "r")
        for expressionData in expressionLevelData.readlines():
            fields = expressionData.strip().split('\t')
            if featureGroupId is not None:
                if fields[3] == featureGroupId:
                    if (expressionLevelId is None or fields[0] ==
                            expressionLevelId):
                        yield self.convertExpressionLevel(fields)
            elif expressionLevelId is None or fields[0] == expressionLevelId:
                yield self.convertExpressionLevel(fields)

# TODO: this is just a first pass stub to get working
# - need to formalize input data
    def convertFeatureGroup(self, record):
        featureGroup = protocol.FeatureGroup()
        featureGroup.id = record[3]
        featureGroup.analysisId = self._rnaQuantificationId
        featureGroup.name = record[3]

        return featureGroup

    def getFeatureGroup(self, featureGroupId):
        """
        for now the feature group data is autogenerated by examining the
        relevant expression data file
        """
        expressionLevelData = open(self._expressionLevelFile, "r")
        for expressionData in expressionLevelData.readlines():
            fields = expressionData.strip().split('\t')
            if fields[3] == featureGroupId:
                yield self.convertFeatureGroup(fields)

    def toProtocolElement(self):
        """
        Converts this rnaQuant into its GA4GH protocol equivalent.
        """
        rnaQuantIterator = self.getRnaQuantification(
            self._rnaQuantificationId)
        rnaQuantData = next(rnaQuantIterator, None)

        protocolElement = protocol.RnaQuantification()
        protocolElement.annotationIds = rnaQuantData.annotationIds
        protocolElement.description = rnaQuantData.description
        protocolElement.id = rnaQuantData.id
        protocolElement.name = rnaQuantData.name
        protocolElement.readGroupId = rnaQuantData.readGroupId

        return protocolElement


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
        analysisId, complexity, exonicFraction, fractionMapped,
        intergenicFraction, intronicFraction
        """
        characterizationData = open(self._characterizationFile, "r")
        quantCharacterization = characterizationData.readline()
        fields = quantCharacterization.split('/t')
        if rnaQuantificationId is None or fields[0] == rnaQuantificationId:
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
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount,
        uniqueSpliceCount
        """
        readCountData = open(self._readCountFile, "r")
        countData = readCountData.readline()
        fields = countData.split('/t')
        if rnaQuantificationId is None or fields[0] == rnaQuantificationId:
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

    def generateExpressionLevel(self):
        """
            Currently just returns default values.
        """
        expressionLevel = protocol.ExpressionLevel

        return expressionLevel

    def getExpressionLevel(self, expressionLevelId, featureGroupId):
        """
        input is tab file with no header.  Columns are:
        annotationId, expression, featureGroupId, id,
        isNormalized, rawReadCount, score, units

        expressionLevelId is not None: return only the specific expressionLevel
        object
        featureGroupId is not None: return all in that group
        """
        expressionLevelData = open(self._expressionLevelFile, "r")
        for expressionData in expressionLevelData.readlines():
            fields = expressionData.strip().split('\t')
            if featureGroupId is not None:
                if fields[3] == featureGroupId:
                    if (expressionLevelId is None or
                            fields[0] == expressionLevelId):
                        yield self.generateExpressionLevel(fields)
            elif expressionLevelId is None or fields[0] == expressionLevelId:
                yield self.generateExpressionLevel(fields)

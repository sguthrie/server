"""
Module responsible for translating feature expression data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import ga4gh.protocol as protocol


class RNASeqResult(object):
    """
    Class representing a single ExpressionAnalysis in the GA4GH data model.
    """
    def __init__(self, expressionAnalysisId, expressionDataPath):
        self._expressionAnalysisId = expressionAnalysisId
        self._expressionAnalysisFile = os.path.join(expressionDataPath, "rnaseq.table")
        self._characterizationFile = os.path.join(expressionDataPath, "dist.table")
        self._readCountFile = os.path.join(expressionDataPath, "counts.table")

    def convertCharacterization(self, record):
        readCharacterization = protocol.Characterization
        readCharacterization.analysisId = record[0]
        readCharacterization.complexity = float(record[1])
        readCharacterization.exonicFraction = float(record[2])
        readCharacterization.fractionMapped = float(record[3])
        readCharacterization.intergenicFraction = float(record[4])
        readCharacterization.intronicFraction = float(record[5])

        return readCharacterization

    def getCharacterization(self, expressionAnalsisId):
        """
        input is tab file with no header.  Columns are:
        analysisId, complexity, exonicFraction, fractionMapped, intergenicFraction, intronicFraction
        """
        characterizationData = open(self._characterizationFile, "r")
        for record in characterizationData.readlines():
            fields = record.split('/t')
            if expressionAnalysisID is None or fields[0] == expressionAnalysisId:
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

    def getReadCounts(self, expressionAnalsisId):
        """
        input is tab file with no header.  Columns are:
        analysisId, multiCount, multiSpliceCount, totalReadCount, uniqueCount, uniqueSpliceCount
        """
        readCountData = open(self._readCountFile, "r")
        for record in readCountData.readlines():
            fields = record.split('/t')
            if expressionAnalysisID is None or fields[0] == expressionAnalysisId:
                yield self.convertReadCounts(fields)

    def convertExpressionAnalysis(self, record):
        expressionAnalysis = protocol.ExpressionAnalysis
        expressionAnalysis.id = record[0]
        expressionAnalysis.annotationIds = record[1].split(',')
        expressionAnalysis.description = record[2]
        expressionAnalysis.name = record[3]
        expressionAnalysis.readGroupId = record[4]

        return expressionAnalysis

    def getExpressionAnalysis(self, expressionAnalysisId):
        """
        input is tab file with no header.  Columns are:
        Id, annotations, description, name, readGroupId
        where annotation is a comma separated list
        """
        expressionAnalysisData = open(self._expressionAnalysisFile, "r")
        for record in expressionAnalysisData.readlines():
            fields = record.strip().split('\t')
            if expressionAnalysisId is None or fields[0] == expressionAnalysisId:
                yield self.convertExpressionAnalysis(fields)

"""
Module responsible for translating reference sequence data into GA4GH native
objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob

import ga4gh.protocol as protocol
import ga4gh.sidegraph as sidegraph


def _makeLimits(start=0, end=None):
    limits = None
    if type(end) is int:
        if start >= 0 and end > start:
            limits = (start, end)
    return limits


def _makeSegment(segmentDict):
    """
    returns a properly formatted protocol.Segment object.
    The input, segmentDict, is a dictionary with keys:
    {sequenceID, start, length, strandIsForward}
    """
    ret = protocol.Segment()
    ret.length = segmentDict["length"]
    ret.start = protocol.Side()
    if segmentDict["strandIsForward"] == sidegraph.SIDEGRAPH_TRUE:
        ret.start.strand = protocol.Strand.POS_STRAND
    else:
        ret.start.strand = protocol.Strand.NEG_STRAND
    ret.start.base = protocol.Position()
    ret.start.base.sequenceId = segmentDict["sequenceID"]
    ret.start.base.position = segmentDict["start"]
    return ret


def _makeJoin(joinDict):
    """
    returns a properly formatted protocol.Join object.
    The input, joinDict, is a dictionary with keys:
    {side1SequenceID, side1Position, side1StrandIsForward,
     side2SequenceID, side2Position, side2StrandIsForward}
    """
    join = protocol.Join()
    # deep pile of objecty dodo
    join.side1 = protocol.Side()
    join.side1.base = protocol.Position()
    if joinDict['side1StrandIsForward'] == sidegraph.SIDEGRAPH_TRUE:
        join.side1.strand = protocol.Strand.POS_STRAND
    else:
        join.side1.strand = protocol.Strand.NEG_STRAND
    join.side1.base.position = joinDict['side1Position']
    join.side1.base.sequenceId = joinDict['side1SequenceID']

    join.side2 = protocol.Side()
    join.side2.base = protocol.Position()
    if joinDict['side2StrandIsForward'] == sidegraph.SIDEGRAPH_TRUE:
        join.side2.strand = protocol.Strand.POS_STRAND
    else:
        join.side2.strand = protocol.Strand.NEG_STRAND
    join.side2.base.position = joinDict['side2Position']
    join.side2.base.sequenceId = joinDict['side2SequenceID']
    return join


class GraphDatabase(object):
    """
    Implements graph reference sets: Reference sequences together with
    joins, as well as some helper methods to deal with complex queries.
    Makes assumption that dataDir used for initialization contains one
    graph topology Sqlite database (.db file) together with any named
    FASTA files containing the actual sequence string data.
    """

    def __init__(self, dataDir):
        self._dataDir = dataDir
        database_files = glob.glob(os.path.join(self._dataDir, "*.sqlite3"))
        assert len(database_files) == 1, "Expects exactly one graph topology Sqlite database (.sqlite3 file) in the data directory: %s" % (self._dataDir)
        self._dbFile = database_files[0]

    def searchReferences(self, referenceSetId=None, start=0, end=None):
        """
        """
        print("searchReferences - referencesGenerator called")
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchReferencesCount()
            referenceDicts = sg.searchReferences(limits)
        references = []
        for rdict in referenceDicts:
            reference = protocol.Reference()
            reference.id = rdict['ID']
            reference.sequenceId = rdict['sequenceID']
            reference.md5checksum = rdict['md5checksum']
            reference.name = rdict['name']
            reference.start = rdict['start']
            reference.length = rdict['length']
            references.append(reference)
        return count, references

    def searchReferenceSets(self, md5checksums=None, accessions=None,
                            assemblyId=None, start=0, end=None):
        """
        Returns a list of dictionaries holding reference set info.
        """
        print("searchReferenceSets - referenceSetsGenerator called")
        # TODO: For now, just returns all reference sets, ignores arguments.
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchReferenceSetsCount()
            referenceSetDicts = sg.searchReferenceSets(limits)
        referenceSets = []
        for rsdict in referenceSetDicts:
            referenceSet = protocol.ReferenceSet()
            referenceSet.id = rsdict['ID']
            referenceSet.includedReferenceSets = []
            referenceSets.append(referenceSet)
        return count, referenceSets

    def searchVariantSets(self, datasetIds=None, start=0, end=None):
        """
        Returns a list of dictionaries holding variant set info.
        """
        print("searchVariantSets - variantSetsGenerator called")
        # TODO: For now, just returns all variant sets, ignores arguments.
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchVariantSetsCount()
            variantSetDicts = sg.searchVariantSets(limits)
        variantSets = []
        for vsdict in variantSetDicts:
            print(vsdict)
            variantSet = protocol.VariantSet()
            variantSet.id = vsdict['ID']
            variantSet.datasetId = ""  # TODO
            variantSet.referenceSetId = ""
            variantSet.metadata = []
            variantSets.append(variantSet)
        return count, variantSets

    def searchVariants(self, start=0, end=None):
        """
        Returns a list of dictionaries holding variant info.
        self.alleleIds = None
        self.alternateBases = None
        self.calls = None
        self.created = None
        self.end = None
        self.id = None
        self.info = {}
        self.names = []
        self.referenceBases = None
        self.referenceName = None
        self.start = None
        self.updated = None
        self.variantSetId = None
        """
        print("searchVariantSets - variantSetsGenerator called")
        # TODO: For now, just returns all variants, since I'm not sure which arguments it would take.
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchVariantsCount()
            variantDicts = sg.searchVariants(limits)
        variants = []
        for vdict in variantDicts:
            variant = protocol.Variant()
            print(vdict)
            variants.append(variant)
        return count, variants

    def searchCallSets(self, datasetIds=None, start=0, end=None):
        """
        """
        print("searchCallSets - callSetsGenerator called")
        # TODO: For now, just returns all variant sets, ignores arguments.
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchCallSetsCount()
            callSetDicts = sg.searchCallSets(limits)
        callSets = []
        for csdict in callSetDicts:
            callSet = protocol.CallSet()
            callSet.id = csdict['ID']
            callSet.sampleId = csdict['sampleID']
            callSet.variantSetIds = csdict['variantSetIds']
            callSet.metadata = []
            callSets.append(callSet)
        return count, callSets

    def searchAlleleCalls(self, alleleId=None,
                          callSetId=None, variantSetId=None,
                          start=0, end=None):
        print("searchAlleleCalls - alleleCallsGenerator called")
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchAlleleCallsCount(alleleId=alleleId,
                                              callSetId=callSetId,
                                              variantSetId=variantSetId)
            alleleCallsDicts = sg.searchAlleleCalls(limits=limits,
                                                    alleleId=alleleId,
                                                    callSetId=callSetId,
                                                    variantSetId=variantSetId)
        alleleCalls = []
        for ac in alleleCallsDicts:
            alleleCall = protocol.AlleleCall()
            alleleCall.callSetId = ac['callSetID']
            alleleCall.alleleId = ac['alleleID']
            alleleCall.totalCopies = ac['ploidy']
            alleleCalls.append(alleleCall)
        return count, alleleCalls

    def searchAlleles(self, variantSetId=None, start=0, end=None):
        print("searchAlleles - allelesGenerator called")
        limits = _makeLimits(start, end)
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchAllelesCount(variantSetId=variantSetId)
            allelesDicts = sg.searchAlleles(limits=limits,
                                            variantSetId=variantSetId)
        alleles = []
        for a in allelesDicts:
            allele = protocol.Allele()
            allele.id = a['ID']
            allele.name = a['name']
            allele.variantSetId = a['variantSetID']
            alleles.append(allele)
        return count, alleles

    def searchSequences(self, referenceSetId=None, variantSetId=None,
                        listBases=False, start=0, end=None):
        """
        Returns a pair (count, sequences) where count
        is the total number of found objects for the search (without limits)
        and sequences is the requested array of protocol.Sequence objects.
        start and end are optional integer arguments that must satisfy
        start < end.

        This method, like its SQL backed siblings,
        does not follow iterator convention.
        """
        print("searchSequences - sequencesGenerator called")
        # TODO search by referenceSetId and variantSetId not yet implemented.
        limits = _makeLimits(start, end)
        count = 0
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchSequencesCount()
            rawSequences = sg.searchSequences(limits)
            sequences = []
            for rawSeq in rawSequences:
                seq = protocol.Sequence()
                seq.id = rawSeq['ID']
                seq.length = rawSeq['length']
                if listBases:
                    seq.bases = sg.getSequenceBases(seq.id)
                sequences.append(seq)
        return count, sequences

    def searchJoins(self, referenceSetId=None, variantSetId=None,
                    sequenceId=None, start=0, end=None):
        """
        Returns a pair (count, joins) where count
        is the total number of found objects for the search (without limits)
        and joins is the requested array of protocol.Join objects.
        start and end are optional integer arguments that must satisfy
        start < end.

        This method, like its SQL backed siblings,
        does not follow iterator convention.
        """
        print("searchJoins - joinsGenerator called")
        # TODO search by referenceSetId and variantSetId not yet implemented.
        # not to mention by sequenceId!
        limits = _makeLimits(start, end)
        count = 0
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            count = sg.searchJoinsCount()
            rawJoins = sg.searchJoins(limits)
        joins = []
        for rj in rawJoins:
            join = _makeJoin(rj)
            joins.append(join)
        return count, joins

    def getSequenceBases(self, sequenceId, start=0, end=None):
        ret = protocol.GetSequenceBasesResponse()
        ret.offset = start
        ret.nextPageToken = None
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            ret.sequence = sg.getSequenceBases(sequenceId, start, end)
        return ret

    def getAllele(self, alleleId):
        """
        Returns a GA4GH formatted allele with the requested ID.
        TODO: Crashes miserably if ID is not valid.
        """
        ret = protocol.Allele()
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            ret.id = alleleId
            alleleData = sg.getAllele(alleleId)
            ret.variantSetId = alleleData['variantSetId']
            ret.name = alleleData['name']
            ret.path = protocol.Path()
            ret.path.segments = map(_makeSegment,
                                    sg.getAllelePathItems(alleleId))
        return ret

    def extractSubgraph(self, seedSequenceId, seedPosition, radius,
                        referenceSetId=None, variantSetId=None):
        """
        Takes a starting (seed) position on a sequence and a radius to
        define a subgraph of all bases and joins reachable by walking
        radius bases (over all allowed joins) from the seed position.

        Returns a pair of arrays: the first of GA4GH-formatted segments,
        the second of GA4GH-formatted joins, which together comprise
        the requested subgraph.
        """
        ret = protocol.ExtractSubgraphResponse()
        with sidegraph.SideGraph(self._dbFile, self._dataDir) as sg:
            segments, joins = sg.getSubgraph(
                seedSequenceId, seedPosition, radius, referenceSetId)
            ret.segments = map(_makeSegment, segments)
            ret.joins = map(_makeJoin, joins)
        return ret

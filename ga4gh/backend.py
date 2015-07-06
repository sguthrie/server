"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import json
import random

import ga4gh.protocol as protocol
import ga4gh.datamodel.references as references
import ga4gh.exceptions as exceptions
import ga4gh.datamodel as datamodel
import ga4gh.datamodel.datasets as datasets


def _parsePageToken(pageToken, numValues):
    """
    Parses the specified pageToken and returns a list of the specified
    number of values. Page tokens are assumed to consist of a fixed
    number of integers seperated by colons. If the page token does
    not conform to this specification, raise a InvalidPageToken
    exception.
    """
    tokens = pageToken.split(":")
    if len(tokens) != numValues:
        msg = "Invalid number of values in page token"
        raise exceptions.BadPageTokenException(msg)
    try:
        values = map(int, tokens)
    except ValueError:
        msg = "Malformed integers in page token"
        raise exceptions.BadPageTokenException(msg)
    return values


def _getVariantSet(request, variantSetIdMap):
    if len(request.variantSetIds) != 1:
        if len(request.variantSetIds) == 0:
            msg = "Variant search requires specifying a variantSet"
        else:
            msg = ("Variant search over multiple variantSets "
                   "not supported")
        raise exceptions.NotImplementedException(msg)
    variantSetId = request.variantSetIds[0]
    try:
        variantSet = variantSetIdMap[variantSetId]
    except KeyError:
        raise exceptions.VariantSetNotFoundException(variantSetId)
    return variantSet


class IntervalIterator(object):
    """
    Implements generator logic for types which accept a start/end
    range to search for the object. Returns an iterator over
    (object, pageToken) pairs. The pageToken is a string which allows
    us to pick up the iteration at any point, and is None for the last
    value in the iterator.
    """
    def __init__(self, request, containerIdMap):
        self._request = request
        self._containerIdMap = containerIdMap
        self._container = self._getContainer()
        self._searchIterator = None
        self._currentObject = None
        self._nextObject = None
        self._searchAnchor = None
        self._distanceFromAnchor = None
        if request.pageToken is None:
            self._initialiseIteration()
        else:
            # Set the search start point and the number of records to skip from
            # the page token.
            searchAnchor, objectsToSkip = _parsePageToken(request.pageToken, 2)
            self._pickUpIteration(searchAnchor, objectsToSkip)

    def _initialiseIteration(self):
        """
        Starts a new iteration.
        """
        self._searchIterator = self._search(
            self._request.start, self._request.end)
        self._currentObject = next(self._searchIterator, None)
        if self._currentObject is not None:
            self._nextObject = next(self._searchIterator, None)
            self._searchAnchor = self._request.start
            self._distanceFromAnchor = 0
            firstObjectStart = self._getStart(self._currentObject)
            if firstObjectStart > self._request.start:
                self._searchAnchor = firstObjectStart

    def _pickUpIteration(self, searchAnchor, objectsToSkip):
        """
        Picks up iteration from a previously provided page token. There are two
        different phases here:
        1) We are iterating over the initial set of intervals in which start
        is < the search start coorindate.
        2) We are iterating over the remaining intervals in which start >= to
        the search start coordinate.
        """
        self._searchAnchor = searchAnchor
        self._distanceFromAnchor = objectsToSkip
        self._searchIterator = self._search(searchAnchor, self._request.end)
        obj = next(self._searchIterator)
        if searchAnchor == self._request.start:
            # This is the initial set of intervals, we just skip forward
            # objectsToSkip positions
            for _ in range(objectsToSkip):
                obj = next(self._searchIterator)
        else:
            # Now, we are past this initial set of intervals.
            # First, we need to skip forward over the intervals where
            # start < searchAnchor, as we've seen these already.
            while self._getStart(obj) < searchAnchor:
                obj = next(self._searchIterator)
            # Now, we skip over objectsToSkip objects such that
            # start == searchAnchor
            for _ in range(objectsToSkip):
                assert self._getStart(obj) == searchAnchor
                obj = next(self._searchIterator)
        self._currentObject = obj
        self._nextObject = next(self._searchIterator, None)

    def next(self):
        """
        Returns the next (object, nextPageToken) pair.
        """
        if self._currentObject is None:
            raise StopIteration()
        nextPageToken = None
        if self._nextObject is not None:
            start = self._getStart(self._nextObject)
            # If start > the search anchor, move the search anchor. Otherwise,
            # increment the distance from the anchor.
            if start > self._searchAnchor:
                self._searchAnchor = start
                self._distanceFromAnchor = 0
            else:
                self._distanceFromAnchor += 1
            nextPageToken = "{}:{}".format(
                self._searchAnchor, self._distanceFromAnchor)
        ret = self._currentObject, nextPageToken
        self._currentObject = self._nextObject
        self._nextObject = next(self._searchIterator, None)
        return ret

    def __iter__(self):
        return self


class ReadsIntervalIterator(IntervalIterator):
    """
    An interval iterator for reads
    """
    def _getContainer(self):
        if len(self._request.readGroupIds) != 1:
            if len(self._request.readGroupIds) == 0:
                msg = "Read search requires a readGroup to be specified"
            else:
                msg = "Read search over multiple readGroups not supported"
            raise exceptions.NotImplementedException(msg)
        readGroupId = self._request.readGroupIds[0]
        try:
            readGroup = self._containerIdMap[self._request.readGroupIds[0]]
        except KeyError:
            raise exceptions.ReadGroupNotFoundException(readGroupId)
        return readGroup

    def _search(self, start, end):
        return self._container.getReadAlignments(
            self._request.referenceId, start, end)

    @classmethod
    def _getStart(cls, readAlignment):
        return readAlignment.alignment.position.position

    @classmethod
    def _getEnd(cls, readAlignment):
        return cls._getStart(readAlignment) + \
            len(readAlignment.alignedSequence)


class VariantsIntervalIterator(IntervalIterator):
    """
    An interval iterator for variants
    """
    def _getContainer(self):
        return _getVariantSet(self._request, self._containerIdMap)

    def _search(self, start, end):
        return self._container.getVariants(
            self._request.referenceName, start, end, self._request.variantName,
            self._request.callSetIds)

    @classmethod
    def _getStart(cls, variant):
        return variant.start

    @classmethod
    def _getEnd(cls, variant):
        return variant.end


class AbstractBackend(object):
    """
    An abstract GA4GH backend.
    This class provides methods for all of the GA4GH protocol end points.
    """
    def __init__(self):
        self._featureIdMap = {}
        self._featureIds = []
        self._referenceSetIdMap = {}
        self._referenceSetIds = []
        self._referenceIdMap = {}
        self._referenceIds = []
        self._requestValidation = False
        self._responseValidation = False
        self._defaultPageSize = 100
        self._maxResponseLength = 2**20  # 1 MiB
        self._datasetIdMap = {}
        self._datasetIds = []

    def _getDatasetFromRequest(self, request):
        if hasattr(request, "datasetIds"):
            if len(request.datasetIds) != 1:
                raise exceptions.NotExactlyOneDatasetException(
                    request.datasetIds)
            datasetId = request.datasetIds[0]
            return self._datasetIdMap[datasetId]
        else:
            # TODO this can go away when all requests have datasetIds
            # instead, we should throw an error here
            datasetId = self._datasetIds[0]
            return self._datasetIdMap[datasetId]

    def getDatasetIds(self):
        """
        Returns a list of datasets in this backend
        """
        return self._datasetIds

    def getDataset(self, datasetId):
        """
        Returns a dataset with id datasetId
        """
        return self._datasetIdMap[datasetId]

    def getReferenceSets(self):
        """
        Returns the list of ReferenceSets in this backend
        """
        return list(self._referenceSetIdMap.values())

    def getReference(self, id_):
        """
        Returns a reference with the given id_
        """
        return self.runGetRequest(self._referenceIdMap, id_)

    def getReferenceSet(self, id_):
        """
        Returns a referenceSet with the given id_
        """
        return self.runGetRequest(self._referenceSetIdMap, id_)

    def getFeature(self, id_):
        """
        Returns a feature with the given id_
        """
        return self.runGetRequest(self._featureIdMap, id_)

    def listReferenceBases(self, id_, requestArgs):
        # parse arguments
        try:
            reference = self._referenceIdMap[id_]
        except KeyError:
            raise exceptions.ObjectWithIdNotFoundException(id_)
        start = 0
        end = datamodel.PysamDatamodelMixin.fastaMax
        if 'start' in requestArgs:
            startString = requestArgs['start']
            try:
                start = int(startString)
            except ValueError:
                raise exceptions.BadRequestIntegerException(
                    'start', startString)
        if 'end' in requestArgs:
            endString = requestArgs['end']
            try:
                end = int(endString)
            except ValueError:
                raise exceptions.BadRequestIntegerException(
                    'end', endString)
        if 'pageToken' in requestArgs:
            pageTokenStr = requestArgs['pageToken']
            start = _parsePageToken(pageTokenStr, 1)[0]
        chunkSize = self._maxResponseLength

        # get reference bases
        gbEnd = min(start + chunkSize, end)
        sequence = reference.getBases(start, gbEnd)

        # determine nextPageToken
        if len(sequence) == chunkSize:
            nextPageToken = start + chunkSize
        elif len(sequence) > chunkSize:
            raise exceptions.ServerError()  # should never happen
        else:
            nextPageToken = None

        # build response
        response = protocol.ListReferenceBasesResponse()
        response.offset = start
        response.sequence = sequence
        response.nextPageToken = nextPageToken
        return response.toJsonString()

    def runGetRequest(self, idMap, id_):
        """
        Runs a get request by indexing into the provided idMap and
        returning a json string of that object
        """
        try:
            obj = idMap[id_]
        except KeyError:
            raise exceptions.ObjectWithIdNotFoundException(id_)
        protocolElement = obj.toProtocolElement()
        jsonString = protocolElement.toJsonString()
        return jsonString

    def runSearchRequest(
            self, requestStr, requestClass, responseClass, objectGenerator):
        """
        Runs the specified request. The request is a string containing
        a JSON representation of an instance of the specified requestClass.
        We return a string representation of an instance of the specified
        responseClass in JSON format. Objects are filled into the page list
        using the specified object generator, which must return
        (object, nextPageToken) pairs, and be able to resume iteration from
        any point using the nextPageToken attribute of the request object.
        """
        self.startProfile()
        try:
            requestDict = json.loads(requestStr)
        except ValueError:
            raise exceptions.InvalidJsonException(requestStr)
        self.validateRequest(requestDict, requestClass)
        request = requestClass.fromJsonDict(requestDict)
        if request.pageSize is None:
            request.pageSize = self._defaultPageSize
        if request.pageSize <= 0:
            raise exceptions.BadPageSizeException(request.pageSize)
        responseBuilder = protocol.SearchResponseBuilder(
            responseClass, request.pageSize, self._maxResponseLength)
        nextPageToken = None
        for obj, nextPageToken in objectGenerator(request):
            responseBuilder.addValue(obj)
            if responseBuilder.isFull():
                break
        responseBuilder.setNextPageToken(nextPageToken)
        responseString = responseBuilder.getJsonString()
        self.validateResponse(responseString, responseClass)
        self.endProfile()
        return responseString

    def searchReadGroupSets(self, request):
        """
        Returns a SearchReadGroupSetsResponse for the specified
        SearchReadGroupSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadGroupSetsRequest,
            protocol.SearchReadGroupSetsResponse,
            self.readGroupSetsGenerator)

    def searchReads(self, request):
        """
        Returns a SearchReadsResponse for the specified
        SearchReadsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadsRequest,
            protocol.SearchReadsResponse,
            self.readsGenerator)

    def searchReferenceSets(self, request):
        """
        Returns a SearchReferenceSetsResponse for the specified
        SearchReferenceSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferenceSetsRequest,
            protocol.SearchReferenceSetsResponse,
            self.referenceSetsGenerator)

    def searchReferences(self, request):
        """
        Returns a SearchReferencesResponse for the specified
        SearchReferencesRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferencesRequest,
            protocol.SearchReferencesResponse,
            self.referencesGenerator)

    def searchVariantSets(self, request):
        """
        Returns a SearchVariantSetsResponse for the specified
        SearchVariantSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantSetsRequest,
            protocol.SearchVariantSetsResponse,
            self.variantSetsGenerator)

    def searchVariants(self, request):
        """
        Returns a SearchVariantsResponse for the specified
        SearchVariantsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantsRequest,
            protocol.SearchVariantsResponse,
            self.variantsGenerator)

    def searchCallSets(self, request):
        """
        Returns a SearchCallSetsResponse for the specified
        SearchCallSetsRequest Object.
        """
        return self.runSearchRequest(
            request, protocol.SearchCallSetsRequest,
            protocol.SearchCallSetsResponse,
            self.callSetsGenerator)

    def searchFeatures(self, request):
        """
        Returns a SearchFeaturesResponse for the specified
        SearchFeaturesRequest object.
        """
        # TODO
        raise exceptions.NotImplementedException("Implement me!")

    def searchRnaQuantification(self, request):
        """
        Returns a SearchRnaQuantificationResponse for the specified
        SearchRnaQuantificationRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchRnaQuantificationRequest,
            protocol.SearchRnaQuantificationResponse,
            self.rnaQuantificationGenerator)

    def searchExpressionLevel(self, request):
        """
        Returns a SearchExpressionLevelResponse for the specified
        SearchExpressionLevelRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchExpressionLevelRequest,
            protocol.SearchExpressionLevelResponse,
            self.expressionLevelGenerator)

    def searchFeatureGroup(self, request):
        """
        Returns a SearchFeatureGroupResponse for the specified
        SearchFeatureGroupRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchFeatureGroupRequest,
            protocol.SearchFeatureGroupResponse,
            self.featureGroupGenerator)

    def searchDatasets(self, request):
        """
        Returns a SearchDatasetsResponse object for the specified
        SearchDatasetsRequest Object.
        """
        return self.runSearchRequest(
            request, protocol.SearchDatasetsRequest,
            protocol.SearchDatasetsResponse,
            self.datasetsGenerator)

    # Iterators over the data hieararchy

    def _topLevelObjectGenerator(self, request, idMap, idList):
        """
        Generalisation of the code to iterate over the objects at the top
        of the data hierarchy.
        """
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = _parsePageToken(request.pageToken, 1)
        while currentIndex < len(idList):
            objectId = idList[currentIndex]
            object_ = idMap[objectId]
            currentIndex += 1
            nextPageToken = None
            if currentIndex < len(idList):
                nextPageToken = str(currentIndex)
            yield object_.toProtocolElement(), nextPageToken

    def datasetsGenerator(self, request):
        """
        Returns a generator over the (dataset, nextPageToken) pairs
        defined by the specified request
        """
        return self._topLevelObjectGenerator(
            request, self._datasetIdMap, self._datasetIds)

    def readGroupSetsGenerator(self, request):
        """
        Returns a generator over the (readGroupSet, nextPageToken) pairs
        defined by the specified request.
        """
        dataset = self._getDatasetFromRequest(request)
        return self._topLevelObjectGenerator(
            request, dataset.getReadGroupSetIdMap(),
            dataset.getReadGroupSetIds())

    def referenceSetsGenerator(self, request):
        """
        Returns a generator over the (referenceSet, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._referenceSetIdMap, self._referenceSetIds)

    def referencesGenerator(self, request):
        """
        Returns a generator over the (reference, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._referenceIdMap, self._referenceIds)

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the specified request.
        """
        dataset = self._getDatasetFromRequest(request)
        return self._topLevelObjectGenerator(
            request, dataset.getVariantSetIdMap(),
            dataset.getVariantSetIds())

    def readsGenerator(self, request):
        """
        Returns a generator over the (read, nextPageToken) pairs defined
        by the specified request
        """
        dataset = self._getDatasetFromRequest(request)
        intervalIterator = ReadsIntervalIterator(
            request, dataset.getReadGroupIdMap())
        return intervalIterator

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined
        by the specified request.
        """
        dataset = self._getDatasetFromRequest(request)
        intervalIterator = VariantsIntervalIterator(
            request, dataset.getVariantSetIdMap())
        return intervalIterator

    def callSetsGenerator(self, request):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined
        by the specified request.
        """
        if request.name is not None:
            raise exceptions.NotImplementedException(
                "Searching over names is not supported")
        dataset = self._getDatasetFromRequest(request)
        variantSet = _getVariantSet(
            request, dataset.getVariantSetIdMap())
        return self._topLevelObjectGenerator(
            request, variantSet.getCallSetIdMap(),
            variantSet.getCallSetIds())

    def rnaQuantificationGenerator(self, request):
        """
        Returns a generator over the (rnaQuantification, nextPageToken) pairs
        defined by the specified request.
        """
        rnaQuantificationId = request.rnaQuantificationId
        if rnaQuantificationId is not None:
            if (rnaQuantificationId in
                    self.getDataset().getRnaQuantificationIds()):
                rnaQuantIds = [rnaQuantificationId]
            else:
                raise exceptions.RnaQuantificationNotFoundException(
                    rnaQuantificationId)
        else:
            rnaQuantIds = self.getDataset().getRnaQuantificationIds()

        return self._topLevelObjectGenerator(
            request,
            self.getDataset().getRnaQuantificationIdMap(), rnaQuantIds)

    def expressionLevelGenerator(self, request):
        expressionLevelId = request.expressionLevelId
        featureGroupId = request.featureGroupId
        rnaQuantificationId = request.rnaQuantificationId
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = _parsePageToken(request.pageToken, 1)
        if rnaQuantificationId is not None:
            rnaQuantificationIds = [rnaQuantificationId]
        else:
            rnaQuantificationIds = self.getDataset().getRnaQuantificationIds()
        rnaQuantifications = self.getDataset().getRnaQuantificationIdMap()
        for rnaQuantId in rnaQuantificationIds:
            rnaQuant = rnaQuantifications[rnaQuantId]
            expressionLevelIterator = rnaQuant.getExpressionLevel(
                expressionLevelId, featureGroupId)
            expressionLevelData = next(expressionLevelIterator, None)
            while (expressionLevelData is not None and
                    currentIndex < self._defaultPageSize):
                nextExpressionLevelData = next(expressionLevelIterator, None)
                nextPageToken = None
                if nextExpressionLevelData is not None:
                    currentIndex += 1
                    nextPageToken = "{}".format(currentIndex)
                expressionLevel = protocol.ExpressionLevel()
                expressionLevel.annotationId = expressionLevelData.annotationId
                expressionLevel.expression = expressionLevelData.expression
                expressionLevel.featureGroupId = (
                    expressionLevelData.featureGroupId)
                expressionLevel.id = expressionLevelData.id
                expressionLevel.isNormalized = expressionLevelData.isNormalized
                expressionLevel.rawReadCount = expressionLevelData.rawReadCount
                expressionLevel.score = expressionLevelData.score
                expressionLevel.units = expressionLevelData.units
                yield expressionLevel, nextPageToken
                expressionLevelData = nextExpressionLevelData

    def featureGroupGenerator(self, request):
        featureGroupId = request.featureGroupId
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = _parsePageToken(request.pageToken, 1)
        rnaQuantifications = self.getDataset().getRnaQuantificationIdMap()
        for rnaQuantId in self.getDataset().getRnaQuantificationIds():
            rnaQuant = rnaQuantifications[rnaQuantId]
            featureGroupIterator = rnaQuant.getFeatureGroup(
                featureGroupId)
            featureGroupData = next(featureGroupIterator, None)
            while (featureGroupData is not None and
                    currentIndex < self._defaultPageSize):
                nextFeatureGroupData = next(featureGroupIterator, None)
                nextPageToken = None
                if nextFeatureGroupData is not None:
                    currentIndex += 1
                    nextPageToken = "{}".format(currentIndex)
                featureGroup = protocol.FeatureGroup()
                featureGroup.analysisId = featureGroupData.analysisId
                featureGroup.created = featureGroupData.created
                featureGroup.description = featureGroupData.description
                featureGroup.id = featureGroupData.id
                featureGroup.info = featureGroupData.info
                featureGroup.name = featureGroupData.name
                featureGroup.updated = featureGroupData.updated
                yield featureGroup, nextPageToken
                featureGroupData = nextFeatureGroupData

    def startProfile(self):
        """
        Profiling hook. Called at the start of the runSearchRequest method
        and allows for detailed profiling of search performance.
        """
        pass

    def endProfile(self):
        """
        Profiling hook. Called at the end of the runSearchRequest method.
        """
        pass

    def validateRequest(self, jsonDict, requestClass):
        """
        Ensures the jsonDict corresponds to a valid instance of requestClass
        Throws an error if the data is invalid
        """
        if self._requestValidation:
            if not requestClass.validate(jsonDict):
                raise exceptions.RequestValidationFailureException(
                    jsonDict, requestClass)

    def validateResponse(self, jsonString, responseClass):
        """
        Ensures the jsonDict corresponds to a valid instance of responseClass
        Throws an error if the data is invalid
        """
        if self._responseValidation:
            jsonDict = json.loads(jsonString)
            if not responseClass.validate(jsonDict):
                raise exceptions.ResponseValidationFailureException(
                    jsonDict, responseClass)

    def setRequestValidation(self, requestValidation):
        """
        Set enabling request validation
        """
        self._requestValidation = requestValidation

    def setResponseValidation(self, responseValidation):
        """
        Set enabling response validation
        """
        self._responseValidation = responseValidation

    def setDefaultPageSize(self, defaultPageSize):
        """
        Sets the default page size for request to the specified value.
        """
        self._defaultPageSize = defaultPageSize

    def setMaxResponseLength(self, maxResponseLength):
        """
        Sets the approximate maximum response length to the specified
        value.
        """
        self._maxResponseLength = maxResponseLength


class EmptyBackend(AbstractBackend):
    """
    A GA4GH backend that contains no data.
    """


class SimulatedBackend(AbstractBackend):
    """
    A GA4GH backend backed by no data; used mostly for testing
    """
    def __init__(self, randomSeed=0, numCalls=1, variantDensity=0.5,
                 numVariantSets=1, numReferenceSets=1,
                 numReferencesPerReferenceSet=1, numAlignments=2):
        super(SimulatedBackend, self).__init__()

        # Datasets
        dataset1 = datasets.SimulatedDataset(
            "simulatedDataset1", randomSeed, numCalls,
            variantDensity, numVariantSets, numAlignments)
        dataset2 = datasets.SimulatedDataset(
            "simulatedDataset2", randomSeed, numCalls,
            variantDensity, numVariantSets, numAlignments)
        self._datasetIdMap[dataset1.getId()] = dataset1
        self._datasetIdMap[dataset2.getId()] = dataset2
        self._datasetIds = sorted(self._datasetIdMap.keys())

        # References
        randomGenerator = random.Random()
        randomGenerator.seed(randomSeed)
        for i in range(numReferenceSets):
            referenceSetId = "referenceSet{}".format(i)
            referenceSetSeed = randomGenerator.getrandbits(32)
            referenceSet = references.SimulatedReferenceSet(
                referenceSetId, referenceSetSeed,
                numReferencesPerReferenceSet)
            self._referenceSetIdMap[referenceSetId] = referenceSet
            for reference in referenceSet.getReferences():
                referenceId = reference.getId()
                self._referenceIdMap[referenceId] = reference
        self._referenceSetIds = sorted(self._referenceSetIdMap.keys())
        self._referenceIds = sorted(self._referenceIdMap.keys())

        # Features
        self._featureIdMap = {}
        self._featureIds = sorted(self._featureIdMap.keys())


class FileSystemBackend(AbstractBackend):
    """
    A GA4GH backend backed by data on the file system
    """
    def __init__(self, dataDir):
        super(FileSystemBackend, self).__init__()
        self._dataDir = dataDir
        # TODO this code is very ugly and should be regarded as a temporary
        # stop-gap until we deal with iterating over the data tree properly.

        # References
        referencesDirName = "references"
        referenceSetDir = os.path.join(self._dataDir, referencesDirName)
        for referenceSetId in os.listdir(referenceSetDir):
            relativePath = os.path.join(referenceSetDir, referenceSetId)
            if os.path.isdir(relativePath):
                referenceSet = references.HtslibReferenceSet(
                    referenceSetId, relativePath)
                self._referenceSetIdMap[referenceSetId] = referenceSet
                for reference in referenceSet.getReferences():
                    referenceId = reference.getId()
                    self._referenceIdMap[referenceId] = reference
        self._referenceSetIds = sorted(self._referenceSetIdMap.keys())
        self._referenceIds = sorted(self._referenceIdMap.keys())

        # Datasets
        datasetDirs = [
            os.path.join(self._dataDir, directory)
            for directory in os.listdir(self._dataDir)
            if os.path.isdir(os.path.join(self._dataDir, directory)) and
            directory != referencesDirName]
        for datasetDir in datasetDirs:
            dataset = datasets.FileSystemDataset(datasetDir)
            self._datasetIdMap[dataset.getId()] = dataset
        self._datasetIds = sorted(self._datasetIdMap.keys())

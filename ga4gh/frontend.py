"""
The Flask frontend for the GA4GH API.

TODO Document properly.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import datetime

import flask
import flask.ext.cors as cors
import humanize
import pydot
import json

import ga4gh
import ga4gh.backend as backend
import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions

MIMETYPE = "application/json"


app = flask.Flask(__name__)


class Version(object):
    """
    A major/minor/revision version tag
    """
    @classmethod
    def parseString(cls, versionString):
        versions = versionString.strip('vV').split('.')
        return Version(*versions)

    def __init__(self, major, minor, revision):
        self.version = (major, minor, revision)

    def __cmp__(self, other):
        return cmp(self.version, other.version)

    def __hash__(self):
        return hash(self.version)

    def __eq__(self, other):
        return self.version == other.version

    def __ne__(self, other):
        return not self.__eq__(other)

    @classmethod
    def getVersionForUrl(cls, versionString):
        """
        Returns the specfied version string in a form suitable for using
        within a URL. This involved prefixing with 'v'.
        """
        ret = versionString
        if not ret.startswith("v"):
            ret = "v{}".format(versionString)
        return ret


class ServerStatus(object):
    """
    Generates information about the status of the server for display
    """

    def __init__(self):
        self.startupTime = datetime.datetime.now()

    def getConfiguration(self):
        """
        Returns a list of configuration (key, value) tuples
        that are useful for users to view on an information page.
        Note that we should be careful here not to leak sensitive
        information. For example, keys and paths of data files should
        not be returned.
        """
        # TODO what other config keys are appropriate to export here?
        keys = [
            'DEBUG', 'REQUEST_VALIDATION', 'RESPONSE_VALIDATION'
        ]
        return [(k, app.config[k]) for k in keys]

    def getPreciseUptime(self):
        """
        Returns the server precisely.
        """
        return self.startupTime.strftime("%H:%M:%S %d %b %Y")

    def getNaturalUptime(self):
        """
        Returns the uptime in a human-readable format.
        """
        return humanize.naturaltime(self.startupTime)

    def getProtocolVersion(self):
        """
        Returns the GA4GH protocol version we support.
        """
        return protocol.version

    def getServerVersion(self):
        """
        Returns the software version of this server.
        """
        return ga4gh.__version__

    def getUrls(self):
        """
        Returns the list of (httpMethod, URL) tuples that this server
        supports.
        """
        urls = []
        rules = app.url_map.iter_rules()
        excluded_methods = ('OPTIONS', 'HEAD')
        excluded_rules = (
            '/', '/flask-api/static/<path:filename>',
            '/static/<path:filename>')
        version = Version.getVersionForUrl(protocol.version)
        for rule in rules:
            for method in rule.methods:
                if (method not in excluded_methods and
                        rule.rule not in excluded_rules):
                    versionedRule = rule.rule.replace(
                        "<version>", version)
                    urls.append((versionedRule, method))
        urls.sort()
        return urls

    def getVariantSets(self):
        """
        Returns the list of variant sets for this server.
        """
        return app.backend.getVariantSets()

    def getReadGroupSets(self):
        """
        Returns the list of ReadGroupSets for this server.
        """
        return app.backend.getReadGroupSets()


def configure(configFile=None, baseConfig="ProductionConfig", extraConfig={}):
    """
    TODO Document this critical function! What does it do? What does
    it assume?
    """
    configStr = 'ga4gh.serverconfig:{0}'.format(baseConfig)
    app.config.from_object(configStr)
    if os.environ.get('GA4GH_CONFIGURATION') is not None:
        app.config.from_envvar('GA4GH_CONFIGURATION')
    if configFile is not None:
        app.config.from_pyfile(configFile)
    app.config.update(extraConfig.items())
    # Setup CORS
    cors.CORS(app, allow_headers='Content-Type')
    app.serverStatus = ServerStatus()
    # Allocate the backend
    # TODO is this a good way to determine what type of backend we should
    # instantiate? We should think carefully about this. The approach of
    # using the special strings __SIMULATED__ and __EMPTY__ seems OK for
    # now, but is certainly not ideal.
    dataSource = app.config["DATA_SOURCE"]
    if dataSource == "__SIMULATED__":
        randomSeed = app.config["SIMULATED_BACKEND_RANDOM_SEED"]
        numCalls = app.config["SIMULATED_BACKEND_NUM_CALLS"]
        variantDensity = app.config["SIMULATED_BACKEND_VARIANT_DENSITY"]
        numVariantSets = app.config["SIMULATED_BACKEND_NUM_VARIANT_SETS"]
        theBackend = backend.SimulatedBackend(
            randomSeed, numCalls, variantDensity, numVariantSets)
    elif dataSource == "__EMPTY__":
        theBackend = backend.EmptyBackend()
    else:
        theBackend = backend.GraphBackend(dataSource)
    theBackend.setRequestValidation(app.config["REQUEST_VALIDATION"])
    theBackend.setResponseValidation(app.config["RESPONSE_VALIDATION"])
    theBackend.setDefaultPageSize(app.config["DEFAULT_PAGE_SIZE"])
    theBackend.setMaxResponseLength(app.config["MAX_RESPONSE_LENGTH"])
    app.backend = theBackend


def getFlaskResponse(responseString, httpStatus=200):
    """
    Returns a Flask response object for the specified data and HTTP status.
    """
    return flask.Response(responseString, status=httpStatus, mimetype=MIMETYPE)


def handleHttpPost(request, endpoint):
    """
    Handles the specified HTTP POST request, which maps to the specified
    protocol handler handpoint and protocol request class.
    """
    if request.mimetype != MIMETYPE:
        raise exceptions.UnsupportedMediaTypeException()
    responseStr = endpoint(request.get_data())
    return getFlaskResponse(responseStr)


def handleHttpOptions():
    """
    Handles the specified HTTP OPTIONS request.
    """
    response = getFlaskResponse("")
    response.headers.add("Access-Control-Request-Methods", "GET,POST,OPTIONS")
    return response


@app.errorhandler(Exception)
def handleException(exception):
    """
    Handles an exception that occurs somewhere in the process of handling
    a request.
    """
    if app.config['DEBUG']:
        app.log_exception(exception)
    serverException = exception
    if not isinstance(exception, exceptions.BaseServerException):
        serverException = exceptions.getServerError(exception)
    responseStr = serverException.toProtocolElement().toJsonString()
    return getFlaskResponse(responseStr, serverException.httpStatus)


def handleFlaskPostRequest(version, flaskRequest, endpoint):
    """
    Handles the specified flask request for one of the POST URLS
    at at the specified version. Invokes the specified endpoint to
    generate a response.
    """
    if Version.parseString(version) != Version.parseString(protocol.version):
        raise exceptions.VersionNotSupportedException()

    if flaskRequest.method == "POST":
        return handleHttpPost(flaskRequest, endpoint)
    elif flaskRequest.method == "OPTIONS":
        return handleHttpOptions()
    else:
        raise exceptions.MethodNotAllowedException()


@app.route('/')
def index():
    return flask.render_template('index.html', info=app.serverStatus)

@app.route('/visual')
def visualize():
    # Get Callsets to find references
    callsets_of_interest = []
    request = protocol.SearchCallSetsRequest()
    request.pageToken = '0'
    while request.pageToken != None:
        responseStr = app.backend.searchCallSets(request.toJsonString())
        response = json.loads(responseStr)
        request.pageToken = response['nextPageToken']
        for callset in response['callSets']:
            if callset['sampleId'].startswith('GRCh'):
                callsets_of_interest.append(callset)

    accepted_callset_ids = [callset['id'] for callset in callsets_of_interest]
    allele_calls = []
    request = protocol.SearchAlleleCallsRequest()
    request.pageToken = '0'
    while request.pageToken != None:
        responseStr = app.backend.searchAlleleCalls(request.toJsonString())
        response = json.loads(responseStr)
        request.pageToken = response['nextPageToken']
        for allele_call in response['alleleCalls']:
            if allele_call['callSetId'] in accepted_callset_ids:
                allele_calls.append(allele_call)


    accepted_allele_ids = [allele_call['alleleId'] for allele_call in allele_calls]
    alleles = []
    request = protocol.SearchAllelesRequest()
    request.pageToken = '0'
    while request.pageToken != None:
        responseStr = app.backend.searchAlleles(request.toJsonString())
        response = json.loads(responseStr)
        request.pageToken = response['nextPageToken']
        for allele in response['alleles']:
            if allele['id'] in accepted_allele_ids:
                alleles.append(allele)

    accepted_allele_paths = []
    sequences_to_find = []
    count, allele_path_items = app.backend.getPathItems()
    for allele_path_item in allele_path_items:
        if allele_path_item['alleleID'] in accepted_allele_ids:
            sequences_to_find.append(allele_path_item['sequenceID'])
            accepted_allele_paths.append(allele_path_item)

    graph = pydot.Dot(graph_type='digraph')
    graph_nodes = {}

    request = protocol.SearchSequencesRequest()
    request.pageToken = '0'
    while request.pageToken != None:
        responseStr = app.backend.searchSequences(request.toJsonString())
        response = json.loads(responseStr)
        request.pageToken = response['nextPageToken']
        for sequence in response['sequences']:
            if sequence['id'] in sequences_to_find:
                name = "ID %s (%s bases)" % (sequence['id'], sequence['length'])
                node = pydot.Node(name)
                graph.add_node(node)
                graph_nodes[sequence['id']] = node

    request = protocol.SearchJoinsRequest()
    request.pageToken = '0'
    while request.pageToken != None:
        responseStr = app.backend.searchJoins(request.toJsonString())
        response = json.loads(responseStr)
        request.pageToken = response['nextPageToken']
        for join in response['joins']:
            n1 = join['side1']['base']['sequenceId']
            n2 = join['side2']['base']['sequenceId']
            if n1 in sequences_to_find and n2 in sequences_to_find:
                node1 = graph_nodes[n1]
                node2 = graph_nodes[n2]
                graph.add_edge(pydot.Edge(node1, node2))

    graph.write_svg('ga4gh/static/tmp_graph.svg')
    return flask.render_template('visual.html', info=app.serverStatus, sequences=response)


@app.route('/<version>/references/<id>', methods=['GET'])
def getReference(version, id):
    raise exceptions.PathNotFoundException()


@app.route('/<version>/references/<id>/bases', methods=['GET'])
def getReferenceBases(version, id):
    raise exceptions.PathNotFoundException()


@app.route('/<version>/referencesets/<id>', methods=['GET'])
def getReferenceSet(version, id):
    raise exceptions.PathNotFoundException()


@app.route('/<version>/callsets/search', methods=['POST'])
def searchCallSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchCallSets)


@app.route('/<version>/allelecalls/search', methods=['POST'])
def searchAlleleCalls(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchAlleleCalls)


@app.route('/<version>/alleles/search', methods=['POST'])
def searchAlleles(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchAlleles)


@app.route('/<version>/readgroupsets/search', methods=['POST'])
def searchReadGroupSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReadGroupSets)


@app.route('/<version>/reads/search', methods=['POST', 'OPTIONS'])
def searchReads(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReads)


@app.route('/<version>/referencesets/search', methods=['POST'])
def searchReferenceSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReferenceSets)


@app.route('/<version>/references/search', methods=['POST'])
def searchReferences(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchReferences)


@app.route('/<version>/variantsets/search', methods=['POST', 'OPTIONS'])
def searchVariantSets(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchVariantSets)


@app.route('/<version>/variants/search', methods=['POST', 'OPTIONS'])
def searchVariants(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchVariants)


@app.route('/<version>/sequences/search', methods=['POST'])
def searchSequences(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchSequences)


@app.route('/<version>/joins/search', methods=['POST'])
def searchJoins(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.searchJoins)


@app.route('/<version>/sequences/<id>/bases', methods=['GET'])
def getSequenceBases(version, id):
    try:
        start = int(flask.request.args.get('start', 0))
    except:
        raise Exception("Start needs to be an integer")
    end = flask.request.args.get('end', None)
    if end is not None:
        try:
            end = int(end)
        except:
            raise Exception("End needs to be an integer")
    responseString = app.backend.getSequenceBases(id, start, end)
    return getFlaskResponse(responseString)


@app.route('/<version>/mode/<mode>', methods=['GET'])
# FIXME Since there's no general GET response mechanism yet,
# this is an utter and unforgivable hack. It may not even be correct,
# as it's unclear to me what a boolean AVRO method should be returning
# (0/1? "True"/"False"? something else entirely?)
def sendsMode(version, mode):
    responseString = 'true' if mode.lower() == 'graph' else 'false'
    return getFlaskResponse(responseString)


@app.route('/<version>/alleles/<id>', methods=['GET'])
def getAllele(version, id):
    responseString = app.backend.getAllele(id)
    return getFlaskResponse(responseString)


@app.route('/<version>/subgraph/extract', methods=['POST'])
def extractSubgraph(version):
    return handleFlaskPostRequest(
        version, flask.request, app.backend.extractSubgraph)


# The below methods ensure that JSON is returned for various errors
# instead of the default, html

@app.errorhandler(404)
def pathNotFoundHandler(errorString):
    return handleException(exceptions.PathNotFoundException())


@app.errorhandler(405)
def methodNotAllowedHandler(errorString):
    return handleException(exceptions.MethodNotAllowedException())

import ga4gh.client as client
import ga4gh.protocol as protocol

#Calls that throw exceptions
#POST 	/v0.6.g/variants/search
#

httpClient = client.HttpClient("http://localhost:8000/v0.6.g", debugLevel=0)
"""
print "==========================================================="
print "CallSets Search"
print "==========================================================="
request = protocol.SearchCallSetsRequest()
response = httpClient.searchCallSets(request)

callsets = []
for thing in response:
    callsets.append(thing)
print len(callsets)
if len(callsets) > 0:
    print callsets[0]

print ""
print "==========================================================="
print "ReferenceSets Search"
print "==========================================================="
request = protocol.SearchReferenceSetsRequest()
response = httpClient.searchReferenceSets(request)

reference_sets = []
for thing in response:
    reference_sets.append(thing)
print len(reference_sets)
if len(reference_sets) > 0:
    print reference_sets[0]

print ""
print "==========================================================="
print "References Search"
print "==========================================================="
request = protocol.SearchReferencesRequest()
response = httpClient.searchReferences(request)

references = []
for thing in response:
    references.append(thing)

print len(references)
if len(references) > 0:
    print references[0]

print ""
print "==========================================================="
print "VariantSets Search"
print "==========================================================="
request = protocol.SearchVariantSetsRequest()
response = httpClient.searchVariantSets(request)

variant_sets = []
for thing in response:
    variant_sets.append(thing)

print len(variant_sets)
if len(variant_sets) > 0:
    print variant_sets[0]

print ""
print "==========================================================="
print "AlleleCalls Search"
print "==========================================================="
request = protocol.SearchAlleleCallsRequest()
response = httpClient.searchAlleleCalls(request)

allele_calls = []
for thing in response:
    allele_calls.append(thing)

print len(allele_calls)
if len(allele_calls) > 0:
    print allele_calls[0]

print ""
print "==========================================================="
print "Alleles Search"
print "==========================================================="
request = protocol.SearchAllelesRequest()
response = httpClient.searchAlleles(request)

alleles = []
for thing in response:
    alleles.append(thing)

print len(alleles)
print alleles[0]
"""
print ""
print "==========================================================="
print "Sequences Search"
print "==========================================================="
request = protocol.SearchSequencesRequest()
response = httpClient.searchSequences(request)

sequences = []
for thing in response:
    sequences.append(thing)

print len(sequences)
print sequences[0]

"""
print ""
print "==========================================================="
print "Joins Search"
print "==========================================================="
request = protocol.SearchJoinsRequest()

response = httpClient.searchJoins(request)

joins = []
for thing in response:
    joins.append(thing)

print len(joins)
print joins[0]


print ""
print "==========================================================="
print "ExtractSubgraph"
print "==========================================================="
request = protocol.ExtractSubgraphRequest()
request.position = protocol.Position()
request.position.sequenceId = int(sequences[0].id)
request.position.position = 1
request.radius = 1

response = httpClient.extractSubgraph(request)

print response
"""

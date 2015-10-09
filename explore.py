import ga4gh.client as client
import ga4gh.protocol as protocol

#Calls that throw exceptions
#POST 	/v0.6.g/variants/search
#

httpClient = client.HttpClient("http://localhost:8000/v0.6.g", debugLevel=0)

"""
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
print "Sequences Search"
print "==========================================================="
request = protocol.SearchSequencesRequest()
response = httpClient.searchSequences(request)

sequences = []
for thing in response:
    sequences.append(thing)

print len(sequences)
print sequences[0]


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

print "==========================================================="
print "CallSets Search - get references"
print "==========================================================="
request = protocol.SearchCallSetsRequest()
response = httpClient.searchCallSets(request)

callsets_of_interest = []
for callset in response:
    if not callset.sampleId.startswith('hu'):
        callsets_of_interest.append(callset)

print ""
print "==========================================================="
print "AlleleCalls Search - get allele calls of references"
print "==========================================================="
request = protocol.SearchAlleleCallsRequest()
response = httpClient.searchAlleleCalls(request)

accepted_callset_ids = [callset.id for callset in callsets_of_interest]
allele_calls = []
for allele_call in response:
    if allele_call.callSetId in accepted_callset_ids:
        allele_calls.append(allele_call)

print ""
print "==========================================================="
print "Alleles Search - get alleles of references"
print "==========================================================="
request = protocol.SearchAllelesRequest()
response = httpClient.searchAlleles(request)

accepted_allele_ids = [allele_call.alleleId for allele_call in allele_calls]
alleles = []
for allele in response:
    if allele.id in accepted_allele_ids:
        print allele
        alleles.append(allele)

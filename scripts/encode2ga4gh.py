#!/usr/bin/env python
from __future__ import division

import os
import sys
import json
import errno
import urllib2
import optparse


class Color:
    """
    Color printing in terminal
    """

    RED = '\033[91m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    NORMAL = '\033[0m'
    YELLOW = '\033[93m'

    @staticmethod
    def red(astr):
        return Color.RED + str(astr) + Color.NORMAL

    @staticmethod
    def green(astr):
        return Color.GREEN + str(astr) + Color.NORMAL

    @staticmethod
    def yellow(astr):
        return Color.YELLOW + str(astr) + Color.NORMAL

    @staticmethod
    def blue(astr):
        return Color.BLUE + str(astr) + Color.NORMAL


def log(message):
    print >>sys.stderr,  "== " + Color.green("{}".format(message))


def get_files_from_host(data, host, output_type, request=False, subset=None):
    file_list = [f for f in data.get('files') if (f.get('output_type')
                 == output_type)]
    if subset:
        file_list = file_list[:subset]
    for f in file_list:
        experiment = f.get('dataset').split('/')[-2]
        replicate = str(f.get('replicate').get('biological_replicate_number'))
        analysisId = "_".join((experiment, replicate))

        if not request:
            yield analysisId
        else:
            href = f.get('href')
            url = "{host}{file}".format(host=host, file=href)
            yield (analysisId, urllib2.urlopen(url))


def get_data_from_host(url, headers, host, output_type, output_folder,
                       description, annotationId, subset=None):
    req = urllib2.Request(url, headers=headers)
    try:
        response = urllib2.urlopen(req)
    except urllib2.URLError as e:
        if hasattr(e, 'reason'):
            print 'We failed to reach a server.'
            print 'Reason: ', e.reason
        elif hasattr(e, 'code'):
            print 'The server couldn\'t fulfill the request.'
            print 'Error code: ', e.code
    else:
        json_data = json.load(response)
        make_dir(output_folder)
        write_rnaseq_tables(get_files_from_host(json_data, host, output_type,
                            subset=subset), description, annotationId,
                            output_folder)
        write_expression_tables(get_files_from_host(json_data, host,
                                output_type, subset=subset, request=True),
                                annotationId, output_folder)
        write_counts_tables(get_files_from_host(json_data, host, output_type,
                            subset=subset, request=True), output_folder)


def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


# TODO: placeholder values need to be calculated then removed
def getCount(expressionId):
    rawCount = 0
    return "%d" % rawCount


# TODO: placeholder values need to be calculated then removed
def getScore(expressionId):
    rawScore = 0.0

    return "%0.2f" % rawScore


def writeRNAQuant(outfile, analysisId, description, annotationId,
                  readGroupId=None):
    """ Using placeholder value for the readGroupId
    """
    if readGroupId is None:
        readGroupId = "readGroupId"

    outline = "\t".join([analysisId, annotationId, description, analysisId,
                         "readGroupId"])
    outfile.write("%s\n" % outline)


def getSamstats(file):
    content = {}
    for line in file.readlines():
        if '|' not in line:
            continue
        key, value = [item.strip() for item in line.split('|')]
        if key == "Number of input reads":
            content["readcount"] = value
        elif key == 'Uniquely mapped reads number':
            content["unique"] = value
        elif key in ['Number of reads mapped to multiple loci',
                     'Number of reads mapped to too many loci']:
            content["multi"] = content.get('multi', 0) + int(value)
        # TODO: STAR does not report unique and multi split-maps, only total
        elif key == 'Number of splices: Total':
            content["msplice"] = content["usplice"] = value

    return content


def writeSamstats(outfile, contents, analysisId):
    outline = "\t".join([analysisId, str(contents["multi"]),
                         contents["msplice"], contents["readcount"],
                         contents["unique"], contents["usplice"]])
    outfile.write("%s\n" % outline)


def getDistribution(file):
    content = {}
    for line in file.readlines():
        split = line.split('\t')
        if split[3] == 'total':
            continue
        content[split[2]] = content.get(split[2], 0) + int(split[3])
    return content


def writeDistribution(outfile, contents, analysisId, fraction):
    outline = "\t".join([analysisId, str(contents["exon"]), fraction,
                         str(contents["intergenic"]), str(contents["intron"])])
    outfile.write("%s\n" % outline)


def write_expression(analysisId, annotationId, quantfile, quantOutfile,
                     tool='RSEM'):
    # RSEM gene expression table header:
    #   gene_id transcript_id(s)    length  effective_length    expected_count
    #   TPM FPKM    pme_expected_count  pme_TPM pme_FPKM    TPM_ci_lower_bound
    #   TPM_ci_upper_bound  FPKM_ci_lower_bound FPKM_ci_upper_bound
    # TODO: placeholder values need to be calculated then removed
    isNormalized = "True"
    units = "TPM"
    # log expression file header
    log(quantfile.readline())
    for expression in quantfile.readlines():
        fields = expression.strip().split("\t")
        expressionLevel = fields[5]
        expressionId = fields[0]
        # TODO: properly handle feature group
        featureGroupId = expressionId
        rawCount = fields[4]
        score = (float(fields[10]) + float(fields[11]))/2
        outline = "\t".join([expressionId, annotationId, expressionLevel,
                             featureGroupId, isNormalized, rawCount,
                             str(score), units])
        quantOutfile.write("%s\n" % outline)


def write_rnaseq_tables(analysisIds, description, annotationId, output_folder):
    log("Writing rnaseq tables")
    for analysisId in analysisIds:
        # create analysis id folder
        make_dir(os.path.join(output_folder, analysisId))

        # output table
        rnaSeqTable = os.path.join(output_folder, analysisId, "rnaseq.table")

        # write rnaseq table
        with open(rnaSeqTable, "w") as rnaQuantFile:
            writeRNAQuant(rnaQuantFile, analysisId, description, annotationId)


def write_counts_tables(data, output_folder):
    log("Writing counts tables")
    for analysisId, samstatsfile in data:
        # output table
        countsTable = os.path.join(output_folder, analysisId, "counts.table")

        samstats = getSamstats(samstatsfile)
        if samstats != {}:
            print(samstats.keys())
            # write mapping stats table
            with open(countsTable, "w") as samOutfile:
                writeSamstats(samOutfile, samstats, analysisId)


def write_dist_tables(data):
    log("Writing distribution tables")
    for analysisId, distfile in data:
        # output table
        distTable = os.path.join(analysisId, "dist.table")

        # write mapping distribution table
        distribution = getDistribution(distfile)
        with open(distTable, "w") as distOutfile:
            writeDistribution(distOutfile, distribution, analysisId,
                              distribution["mapped"])


def write_expression_tables(data, annotationId, output_folder):
    log("Writing gene expression tables")
    for analysisId, quantfile in data:
        # output table
        expTable = os.path.join(output_folder, analysisId, "expression.table")

        # write expression table
        print("processing {}".format(analysisId))
        with open(expTable, "w") as quantOutfile:
            write_expression(analysisId, annotationId, quantfile, quantOutfile)


def makeParser(usage):
    """
        Default values will download and parse 4 quantifications from the
        ENCODE Evaluation dataset as an example set of test data.
    """

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--dataset", dest="dataset")
    parser.add_option("--description", dest="description")
    parser.add_option("--annotation", dest="annotationId")
    parser.add_option("--file_limit", type="int", dest="subset")

    parser.set_defaults(dataset="ENCSR000AJW", annotationId="Gencodev16",
                        description="RNAseq data from ENCODE evaluation",
                        subset=4)

    return parser


def main(argv):

    usage = "Usage: {} <data-folder>".format(argv[0])
    if len(argv) != 2:
        print(Color.red(usage))
        sys.exit(1)

    parser = makeParser(usage)
    (options, args) = parser.parse_args(argv[1:])

    host = "https://www.encodeproject.org"
    dataset = options.dataset
    url = "{host}/datasets/{dataset}/?frame=embedded".format(host=host,
                                                             dataset=dataset)
    headers = {
        'Accept': 'application/json; charset=utf-8'
    }

    data_type = "gene quantifications"
    data_folder = argv[1]
    rna_folder = "rnaQuant"
    output_folder = os.path.join(data_folder, rna_folder)
    subset = options.subset

    log("Downloading GA4GH test dataset - RNA Quantification API")
    print("ENCODE dataset: {}".format(Color.blue(dataset)))
    print("data type:      {}".format(Color.blue(data_type)))
    print("subset size:    {}".format(Color.blue(subset)))
    print("output folder:  {}".format(Color.blue(output_folder)))
    get_data_from_host(url, headers, host, data_type, output_folder,
                       options.description, options.annotationId, subset)

    log("DONE")


if __name__ == '__main__':
    main(sys.argv)

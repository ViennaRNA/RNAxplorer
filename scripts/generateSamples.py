"""
General pipeline for diverse structure sample set generation
using RNAxplorer
"""

import sys, re, RNA


def readFasta(filename):
    """
    Read a FASTA file and extract records as pairs of
    id, and sequence
    """
    res = []
    counter = 0
    is_fasta  = 0
    record_id = "seq_" + `counter`
    sequence = ""
    for l in open(filename):
        l.rstrip('\n')

        # start actual parsing
        if l.startswith(">"):
            match = re.search('>\s*(\S+)', l)
            if match:
                is_fasta  = 1
                record_id = match.group(1)
                sequence  = ""
        elif not l.startswith("#")  and not l.startswith(";"):
            match = re.search('([ACGUTacgutnN]+)', l)
            if match:
                sequence += match.group(1)

        if not is_fasta:
            res.append((record_id, sequence))
            sequence = ""
            counter += 1
            record_id = "seq_" + `counter`

    return res

def callRNAxplorer(seq, ref_struct1, ref_struct2, n=100):
    """
    Call RNAxplorer using ref_struct1, ref_struct2 as
    structures that are made equally probable within
    the ensemble. Tell RNAxplorer to generate n samples

    At this point, we might have redundancy in the sample
    set which needs to be removed later
    """
    structures = []

    print "Generating %d samples with RNAxplorer using\n%s\n%s\n%s" % (n, seq, ref_struct1, ref_struct2)

    return structures


def addToClusters(clusters, structures):
    """
    Group the provided structures into
    clusters, and for each cluster compute
    a representative, e.g. the centroid, or
    MFE

    We can either utilize the hierarchical clustering method
    that we employed before, or, as Maria suggested, create
    gradient basins at this stage.
    """
    new_clusters = clusters

    return new_clusters


def mergeStructures(samples, new_structures):
    """
    Merge new_structures into a pre-existing
    sample set non-redundantly
    """

    # to be implemented
    samples.append(new_structures)

    return samples


def generateSamples(seq, mfe_struct, mfe, iterations=0):
    """
    Generate a diverse structure sample set
    by iteratively calling RNAxplorer
    """

    # Initial step:
    # start with mfe_struct and open chain as references
    # to retrieve first sample set
    openchain = "." * len(seq)
    samples = callRNAxplorer(seq, mfe_struct, openchain)

    # add new clusters from current structure sample set
    clusters = []
    clusters = addToClusters(clusters, samples)

    # here comes the iterative step
    for i in range(iterations):
        interesting_clusters = getInterestingClusters(clusters)
        for representative in intersting_clusters:
            # create more samples by making mfe_struct and representative
            # equally probable
            new_samples = callRNAxplorer(seq, mfe_struct, representative)
            clusters = addToClusters(clusters, new_samples)

    return clusters


if __name__ == "__main__":
    inputfile = sys.argv[1]
    records = readFasta(inputfile)

    for fasta_id, seq in records:
        print "...processing new Record:\n\nID: %s\n%s" % (fasta_id, seq)

        # compute MFE and corresponding structure
        mfe_struct, mfe = RNA.fold(seq)
        print "%s (%6.2f)" % (mfe_struct, mfe)

        # generate diverse sample set
        clusters = generateSamples(seq,mfe_struct,mfe)


        print "\n...done\n"

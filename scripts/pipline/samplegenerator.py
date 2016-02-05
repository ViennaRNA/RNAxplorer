"""
General pipeline for diverse structure sample set generation
using RNAxplorer
"""

import sys, re, RNA, subprocess, operator
from clusteralgorithms.two_d_flooder import WatershedFlooder

class Matrix2D:
    """
    contains all the information about an rna_xplorer call.
    """
    def __init__(self, sequence, samples, ss1, ss2):
        """
        Args:
            sequence: rna sequence [ACGTacgt]
            samples: list with entries of the form [bpdist to ref1, pbdist to ref2, free energy, secondary strucutre]
            ss1: first reference structure
            ss2: second reference structure
        """
        self.seq = sequence
        self.ss1 = ss1
        self.ss2 = ss2
        self.samples = samples

def callRNAxplorer(seq, ref_struct1, ref_struct2, n=500):
    """
    Call RNAxplorer using ref_struct1, ref_struct2 as
    structures that are made equally probable within
    the ensemble. Tell RNAxplorer to generate n samples

    returns a Matrix2D object which contains all unique samples and the inputdata of the xplorer call.
    """
    xplorerSamples = []
    unique = " | sort -k5 | uniq"
    RNAxplorer = "RNAxplorer -M SM -e \"MSF\" --betaScale=1.2 -i " + str(n) + unique

    sequenceAndStructures = seq + "\n" + ref_struct1 + "\n" + ref_struct2
    # Run RNAxplorer
    to_call = "echo -e \"" + sequenceAndStructures + "\" | " + RNAxplorer
    result = subprocess.check_output(to_call, shell=True)
    result = result.splitlines()
    for line in result:   
        match = re.search('^(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s+\(\d+\)\s+([\(\)\.]+)$', line)
        if match:
            k = match.group(1)  # distance to reference structure one
            l = match.group(2)  # distance to reference structure two
            e = match.group(3)  # free energy of the structure
            s = match.group(4)  # the structure
            entry = [int(k), int(l), float(e), s]
            xplorerSamples.append(entry)
    
    print "Generating %d samples with RNAxplorer using\n%s\n%s\n%s" % (n, seq, ref_struct1, ref_struct2)
    
    x = Matrix2D(seq, xplorerSamples, ref_struct1, ref_struct2)
    return x


def addToClusters(clusters, structures):
    """
    Group the provided structures into clusters,
    and for each cluster compute a representative,
    e.g. the centroid, or MFE, and possibly its
    ensemble diversity and partition function,
    i.e. whatever property we might be interested in

    We can either utilize the hierarchical clustering
    method that we employed before, or, as Maria
    suggested, create gradient basins at this stage.
    """
    new_clusters = clusters

    # this loop is just a placeholder
    # we actually need to assign each structure to
    # a corresponding cluster instead
    for s in structures:
        new_clusters.append(s)

    return new_clusters


def mfeStructure(seq, cluster):
    mfe = sys.float_info.max
    mfe_s = ""
    for s, i in cluster:
        e = RNA.energy_of_struct(seq, s)
        if e < mfe:
            mfe = e
            mfe_s = s
    return mfe_s

def getInterestingClusters(seq, clusters):
    """
    Extract a list of interesting structures that
    represent the provided clusters. These structures
    are then used in the next sampling iteration to
    emphasize a particular region in the landscape.

    This could be a cluster with low sample number,
    or one that has the highest distance to all the
    others. Open for suggestions here...
    """
    
    """
    #lets try all mfe representatives of the clusters.
    selection = []
    
    for i, b in enumerate(clusters):
        minE = sys.maxint
        minS = ""
        for (s,idx) in b:
            energy = RNA.energy_of_struct(seq, s)
            if energy < minE:
                minS = s
                minE = energy
        selection.append(minS)
    """
    
    # try two smallest basins
    selection = []
    for c in clusters: 
        length = len(c)
        selection.append((length, c))
    selection.sort(key=operator.itemgetter(0))
    
    
    b1 = selection[0][1]
    b2 = selection[1][1]
    s1 = mfeStructure(seq, b1)
    s2 = mfeStructure(seq, b2)
    
    return [s1, s2]


def get2DBasins(samples, sequence, ref1, ref2):
    """
    Given two reference structures ref1, and ref2, generate
    the 2D projection of the provided sample set of structures.
    Using the MFE representatives of each distance class,
    determine the basins of attraction in this projection,
    and associate each structure to such basin.
    Basin detection might be implemented using a flooding
    technique that stops as soon as a saddle point in the
    projection is encountered
    
    This function returns a list of lists of (s,idx) pairs, where
    s is a secondary structure in dot-bracket notation, and idx
    is the corresponding index in the samples list
    """

    return WatershedFlooder.doflooding(samples)

def uniq_list(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def generateSamples(seq, mfe_struct, reference_stack):
    """
    Generate a diverse structure sample set by calling
    RNAxplorer for each pair of reference structures
    that is supplied

    After calling RNAxplorer we extract 2D basins from
    the projection and retrieve a set of new interesting
    reference structures from them.
    """

    # The global storage of all structures sampled in this run.
    sample_set = set()

    # stack to hold interesting new reference structures
    new_reference_stack = []


    # now process all the structures in the reference_stack
    # to obtain more structure samples. 
    while len(reference_stack) >= 2:
        # at this point, pop two structures at a time to get
        # new references and only default to using mfe_struct
        # as first reference if we have just a single new_ref left
        new_ref1 = reference_stack.pop()
        if len(reference_stack) > 0:
            new_ref2 = reference_stack.pop()
        else:
            new_ref2 = mfe_struct

        # call RNAxplorer with our reference structures
        xplorerData = callRNAxplorer(seq, new_ref1, new_ref2)
        new_samples = xplorerData.samples

        # determine which structures belong to the same region
        # of interest, according to current 2D projection.
        basins_2D = get2DBasins(new_samples, seq, new_ref1, new_ref2)

        # cluster the structures within each basin to
        # obtain potentially interesting new reference
        # structures
        nextClusters = getInterestingClusters(seq, basins_2D)
        for c in nextClusters:
            new_reference_stack.append(c)    
                  
            # basin_clusters = getClusters(b)
            # s = getInterestingClusters(basin_clusters)
            # for c in s:
            #    new_reference_stack.append(c)

        # finally, we add the new_samples to our sample_set
        for s in new_samples:
            structure = {s[3]}
            sample_set.update(structure)

    return sample_set, new_reference_stack


def mainloop(seq, mfe_struct, mfe, iterations=0):
    """
    This is the mainloop that iteratively calls
    generateSamples() to grow the global_structures list
    """
    # the global storage of all structures sampled so far.
    # this variable is called clusters, however, we actually
    # don't require the structures to be clustered here.
    global_structures = set()

    # Initial step:
    # start with mfe_struct and open chain as references
    # to retrieve first sample set
    ref_structs = []
    openchain = "." * len(seq)
    ref_structs.append(mfe_struct)
    ref_structs.append(openchain)

    while iterations >= 0 and len(ref_structs) > 0:
        # get new samples and potentially new reference structures
        samples, new_ref_structs = generateSamples(seq, mfe_struct, ref_structs)
        
        # add new samples to set of global structures
        global_structures.update(samples)
        # global_structures = addToClusters(global_structures, samples) 
        
        # assign new reference structures for next iteration   
        ref_structs = new_ref_structs
        iterations = iterations - 1

    return global_structures


#/usr/bin/env python2
#

from __future__ import print_function
import sys
import math

import RNA
import RNAxplorer

from pipeline.clusteralgorithms.diana import DIANA

num_iter = 15
num_samples = 1000

# this is SV11
sequence = "GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA"

"""
Error print
"""
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def RNAlocmin_output(seq, m):
    print("     %s" % seq)

    for i,s in enumerate(sorted(m.keys(), key=lambda x: m[x]['energy'])):
        print("%4d %s %6.2f %6d" % (i, s, m[s]['energy'], m[s]['count']))


def mfeStructure(seq, cluster):
    mfe = sys.float_info.max
    mfe_s = ""
    for s in cluster:
        e = RNA.energy_of_struct(seq, s)
        if e < mfe:
            mfe = e
            mfe_s = s
    return mfe_s


def selectRepresentatives(seq, clusters):
    """
    extract mfe structures from clusters
    @param seq - string - the RNA sequence.
    @return list - structures in dot-bracket notation
    """
    representatives = []
    for c in clusters:
        # select mfe or TODO: centroid structure.
        s = mfeStructure(seq, c)
        representatives.append(s)
    return representatives


# init random number generator in RNAlib
RNA.init_rand()

# prepare RNAlib fold_compound
md = RNA.md()
md.uniq_ML = 1
md.compute_bpp = 0

fc = RNA.fold_compound(sequence, md, RNA.OPTION_PF)
fc_base = RNA.fold_compound(sequence, md)

# compute MFE structure
(ss, mfe) = fc.mfe()


minima = dict()
repulsed_structures = dict()

for it in range(0, num_iter):
    eprint("iteration %d" % it)

    # fill partition function DP matrices 
    fc.pf()

    for i in range(0, num_samples):
        # sample structure
        s = fc.pbacktrack()
        # perform gradient walk from sample to determine direct local minimum
        pt = RNA.IntVector(RNA.ptable(s))
        fc_base.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
        ss = RNA.db_from_ptable(list(pt))
        if ss not in minima:
            minima[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
        else:
            minima[ss]['count'] = minima[ss]['count'] + 1

    # find out which local minima we've seen the most
    repell_struct = max(minima.iterkeys(), key=(lambda a: minima[a] if a not in repulsed_structures else 0))
    repell_en = - minima[repell_struct]['energy'] / 3.0

    eprint("repelling the following struct")
    eprint("%s (%6.2f)" % (repell_struct, repell_en))

    repulsed_structures[repell_struct] = 1

    # add repell_struct to repulsion potential
    # and lift it half teh way up to 0-base line
    if it < num_iter - 1:
        RNAxplorer.add_repulsion(fc, repell_struct, repell_en)


initial_minima = [ s for s in minima.keys() ]

maxDiameterThreshold = 0
maxAverageDiameterThreshold = 6

clusters = DIANA.doClustering(initial_minima, maxDiameterThreshold, maxAverageDiameterThreshold)
#DIANA.printClusters(clusters)

eprint("done clustering %d structures into %d clusters" % (len(initial_minima), len(clusters)))

representatives = selectRepresentatives(sequence, clusters)

# create a minima-like dict for output printing
red_minima = dict()
for r in representatives:
    red_minima[r] = {'count': 1, 'energy' : fc_base.eval_structure(r)}

# produce RNAlocmin - like output
#RNAlocmin_output(sequence, minima)
RNAlocmin_output(sequence, red_minima)


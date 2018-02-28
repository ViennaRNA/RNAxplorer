#/usr/bin/env python2
#

from __future__ import print_function
import sys
import math
import argparse

import RNA
import RNAxplorer

from pipeline.clusteralgorithms.diana import DIANA

num_iter = 150
num_samples = 1000
kt_fact = 1
min_explore_min_percent = 10
exploration_factor = 1.2
do_clustering = True
verbose = False
lmin_file = "local_minima.txt"

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


"""
Start argument parsing
"""
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--verbose", help="Be verbose", action="store_true")
parser.add_argument("-s", "--sequence", type=str, help="Input sequence")
parser.add_argument("-i", "--iterations", type=int, help="Number of iterations")
parser.add_argument("-n", "--num-samples", type=int, help="Number of samples per iteration")
parser.add_argument("-f", "--exploration-factor", type=float, help="Exploration factor")
parser.add_argument("--min-exploration-percent", type=float, help="Minimum exploration percentage before adding new repelled structures")
parser.add_argument("-c", "--cluster", help="Cluster resulting local minima to reduce effective size", action="store_true")
parser.add_argument("--no-cluster", help="Do not cluster resulting local minima to reduce effective size", action="store_true")
parser.add_argument("--lmin-file", type=str, help="Output filename for local minima")


args = parser.parse_args()

if args.verbose:
    verbose = True

if args.sequence:
    sequence = args.sequence

if args.iterations:
    num_iter = args.iterations

if args.num_samples:
    num_samples = args.num_samples

if args.no_cluster:
    do_clustering = False

if args.exploration_factor:
    exploration_factor = args.exploration_factor

if args.min_exploration_percent:
    min_explore_min_percent = args.min_exploration_percent

if args.lmin_file:
    lmin_file = args.lmin_file


"""
Do main stuff
"""
# init random number generator in RNAlib
RNA.init_rand()

# prepare RNAlib fold_compound
md = RNA.md()
md.uniq_ML = 1
md.compute_bpp = 0

kT = RNA.exp_param(md).kT

fc = RNA.fold_compound(sequence, md, RNA.OPTION_PF)
fc_base = RNA.fold_compound(sequence, md)

# compute MFE structure
(ss, mfe) = fc.mfe()


minima = dict()
repulsed_structures = dict()

last_reference_id = -1
last_reference_weight = -1

for it in range(0, num_iter):
    if verbose:
        eprint("iteration %d" % it)

    # fill partition function DP matrices 
    fc.pf()

    new_minima = 0

    for i in range(0, num_samples):
        # sample structure
        s = fc.pbacktrack()
        # perform gradient walk from sample to determine direct local minimum
        pt = RNA.IntVector(RNA.ptable(s))
        fc_base.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
        ss = RNA.db_from_ptable(list(pt))
        if ss not in minima:
            minima[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
            new_minima = new_minima + 1
        else:
            minima[ss]['count'] = minima[ss]['count'] + 1

    if verbose:
        eprint("new minima: %d vs. present minima: %d" % (new_minima, len(minima)))

    if it < num_iter - 1:
        if (last_reference_id != -1) and (new_minima == 0 or (len(minima) / new_minima) > min_explore_min_percent) and (last_reference_weight < -mfe):
            if verbose:
                eprint("increasing strength of last reference (id: %d, from %6.2f to %6.2f)" % (last_reference_id, last_reference_weight, last_reference_weight * exploration_factor))

            RNAxplorer.change_repulsion(fc, last_reference_id, last_reference_weight * exploration_factor)
            last_reference_weight = last_reference_weight * exploration_factor

        else:
            # find out which local minima we've seen the most
            repell_struct = max(minima.iterkeys(), key=(lambda a: minima[a] if a not in repulsed_structures else 0))
            repell_en = kt_fact * kT / 1000.
            last_reference_weight = repell_en

            if verbose:
                eprint("repelling the following struct (last id: %d)" % last_reference_id)
                eprint("%s (%6.2f)" % (repell_struct, repell_en))

            repulsed_structures[repell_struct] = 1

            last_reference_id = RNAxplorer.add_repulsion(fc, repell_struct, repell_en)


# save local minima to file
f = open(lmin_file, 'w')
f.write("     %s\n" % sequence)
for i,s in enumerate(sorted(minima.keys(), key=lambda x: minima[x]['energy'])):
        f.write("%4d %s %6.2f %6d\n" % (i, s, minima[s]['energy'], minima[s]['count']))
f.close()

if do_clustering:
    initial_minima = [ s for s in minima.keys() ]

    maxDiameterThreshold = 0
    maxAverageDiameterThreshold = 6

    clusters = DIANA.doClustering(initial_minima, maxDiameterThreshold, maxAverageDiameterThreshold)
    #DIANA.printClusters(clusters)

    if verbose:
        eprint("done clustering %d structures into %d clusters" % (len(initial_minima), len(clusters)))

    representatives = selectRepresentatives(sequence, clusters)

    # create a minima-like dict for output printing
    red_minima = dict()
    for r in representatives:
        red_minima[r] = {'count': 1, 'energy' : fc_base.eval_structure(r)}

    # produce RNAlocmin - like output
    #RNAlocmin_output(sequence, minima)
    RNAlocmin_output(sequence, red_minima)


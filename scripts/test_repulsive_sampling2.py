#/usr/bin/env python2
#

from __future__ import print_function
import sys
import math
import argparse

import RNA
import RNAxplorer

from pipeline.clusteralgorithms.diana import DIANA

num_iter = 100
num_samples = 10000
kt_fact = 1
min_explore_min_percent = 10
exploration_factor = 1.2
do_clustering = False
fake_2D_file = True
# switch to select penalization of base pairs instead of computing complicated structure penalties
penalize_base_pairs = True
verbose = False
debug   = False
lmin_file = "local_minima.txt"
nonredundant_sample_file = False

# this is SV11
sequence = "GGGCACCCCCCUUCGGGGGGUCACCUCGCGUAGCUAGCUACGCGAGGGUUAAAGCGCCUUUCUCCCUCGCGUAGCUAACCACGCGAGGUGACCCCCCGAAAAGGGGGGUUUCCCA"
sv11_mfe="(((.(((((((((((((((((((((((((((.(.((((((((((((((..((((...))))..)))))))))))))).).)))))))))))))))))))..))))))))..)))."
sv11_meta="(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((........)))))))).((((((((.....))))))))....."

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
parser.add_argument("-d", "--debug", help="Be even more verbose", action="store_true")
parser.add_argument("-s", "--sequence", type=str, help="Input sequence")
parser.add_argument("-i", "--iterations", type=int, help="Number of iterations")
parser.add_argument("-n", "--num-samples", type=int, help="Number of samples per iteration")
parser.add_argument("-f", "--exploration-factor", type=float, help="Exploration factor")
parser.add_argument("--min-exploration-percent", type=float, help="Minimum exploration percentage before adding new repelled structures")
parser.add_argument("-c", "--cluster", help="Cluster resulting local minima to reduce effective size", action="store_true")
parser.add_argument("--lmin-file", type=str, help="Output filename for local minima")
parser.add_argument("--nonred-file", type=str, help="Input filename for nonredundant samples")


args = parser.parse_args()

if args.verbose:
    verbose = True

if args.debug:
    debug = True

if args.sequence:
    sequence = args.sequence

if args.iterations:
    num_iter = args.iterations

if args.num_samples:
    num_samples = args.num_samples

if args.cluster:
    do_clustering = True

if args.exploration_factor:
    exploration_factor = args.exploration_factor

if args.min_exploration_percent:
    min_explore_min_percent = args.min_exploration_percent

if args.lmin_file:
    lmin_file = args.lmin_file

if args.nonred_file:
    nonredundant_sample_file = args.nonred_file



sc_data = {
  'base_pairs': {},
  'weights': {},
}


def store_basepair_sc(data, structure, weight):
    pt = RNA.ptable(structure)
    # count number of pairs in structure to repell
    cnt = 0
    for i in range(1, len(pt)):
        if pt[i] > i:
            cnt = cnt + 1

    if cnt > 0:
        weight = weight / cnt

    # add repulsion
    for i in range(1, len(pt)):
        if pt[i] > i:
            key = (i, pt[i])
            if key not in data['base_pairs']:
                data['base_pairs'][key] = 1
                data['weights'][key] = weight
                if debug:
                    print("adding pair (%d,%d) with weight %g" % (i, pt[i], weight))
            else:
                    data['base_pairs'][key] = data['base_pairs'][key] + 1
                    data['weights'][key] = data['weights'][key] + weight
                    if debug:
                        print("We've seen this pair (%d,%d) before! Increasing its repellent potential to %g)!" % (i, pt[i], data['weights'][key]))


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

fc = RNA.fold_compound(sequence, md)
fc_base = RNA.fold_compound(sequence, md)

# compute MFE structure
(ss, mfe) = fc.mfe()


minima = dict()
repulsed_structures = dict()


for it in range(0, num_iter):
    if verbose:
        eprint("iteration %d" % it)

    #fc = RNA.fold_compound(sequence, md)
    fc.sc_remove()

    for k in sc_data['weights'].keys():
        i = k[0]
        j = k[1]
        fc.sc_add_bp(i, j, sc_data['weights'][k])

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
        # find out which local minima we've seen the most
        repell_struct = max(minima.iterkeys(), key=(lambda a: minima[a] if a not in repulsed_structures else 0))
        repell_en = kt_fact * kT / 1000.

        if verbose:
            eprint("repelling the following struct\n%s (%6.2f)" % (repell_struct, repell_en))

        repulsed_structures[repell_struct] = 1

        store_basepair_sc(sc_data, repell_struct, repell_en)


# save local minima to file
f = open(lmin_file, 'w')
f.write("     %s\n" % sequence)
for i,s in enumerate(sorted(minima.keys(), key=lambda x: minima[x]['energy'])):
        f.write("%4d %s %6.2f %6d\n" % (i, s, minima[s]['energy'], minima[s]['count']))
f.close()

if fake_2D_file:
    distances = [ [ None for j in range(0, 200) ] for i in range(0, 200) ];

    for s in minima.keys():
        d1 = RNA.bp_distance(sv11_mfe, s)
        d2 = RNA.bp_distance(sv11_meta, s)
        if not distances[d1][d2] or minima[s]['energy'] < distances[d1][d2]:
            distances[d1][d2] = minima[s]['energy']

    f = open("sv11_fake.2D.out", 'w')
    f.write("%s\n%s (%6.2f)\n%s\n%s\n\n\n" % (sequence, sv11_mfe, mfe, sv11_mfe, sv11_meta))
    for i in range(0, 200):
        for j in range(0, 200):
            if distances[i][j] != None:
                f.write("%d\t%d\ta\tb\tc\t%6.2f\n" % (i, j, distances[i][j]))
    f.close()


# read a list of sample structures and produce list of local minima for it
if nonredundant_sample_file:
    lmin_nonred_file = "local_minima_nonred.txt"
    nonredundant_samples = []
    with open(nonredundant_sample_file) as f:
        nonredundant_samples = f.readlines()

    nonredundant_samples = [x.strip() for x in nonredundant_samples]

    nonredundant_samples.pop(0)

    nonredundant_minima = dict()

    for s in nonredundant_samples:
        pt = RNA.IntVector(RNA.ptable(s))
        fc_base.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
        ss = RNA.db_from_ptable(list(pt))
        if ss not in nonredundant_minima:
             nonredundant_minima[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
             new_minima = new_minima + 1
        else:
             nonredundant_minima[ss]['count'] = nonredundant_minima[ss]['count'] + 1

    f = open(lmin_nonred_file, 'w')
    f.write("     %s\n" % sequence)
    for i,s in enumerate(sorted(nonredundant_minima.keys(), key=lambda x: nonredundant_minima[x]['energy'])):
        f.write("%4d %s %6.2f %6d\n" % (i, s, nonredundant_minima[s]['energy'], nonredundant_minima[s]['count']))
    f.close()

    if fake_2D_file:
        distances = [ [ None for j in range(0, 200) ] for i in range(0, 200) ];

        for s in nonredundant_minima.keys():
            d1 = RNA.bp_distance(sv11_mfe, s)
            d2 = RNA.bp_distance(sv11_meta, s)
            if not distances[d1][d2] or nonredundant_minima[s]['energy'] < distances[d1][d2]:
                distances[d1][d2] = nonredundant_minima[s]['energy']

        f = open("sv11_fake_nonred.2D.out", 'w')
        f.write("%s\n%s (%6.2f)\n%s\n%s\n\n\n" % (sequence, sv11_mfe, mfe, sv11_mfe, sv11_meta))
        for i in range(0, 200):
            for j in range(0, 200):
                if distances[i][j] != None:
                    f.write("%d\t%d\ta\tb\tc\t%6.2f\n" % (i, j, distances[i][j]))
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


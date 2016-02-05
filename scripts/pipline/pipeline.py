#!/usr/bin/env python

import RNA
import samplegenerator
from clusteralgorithms.diana import DIANA
from rates import Rates
import filewriter
import os
import operator

def selectRepresentatives(clusters):
    representatives = []
    for c in clusters:
        # TODO: select mfe or centroid structure.
        s = c[0]
        representatives.append(s)
    return representatives
    

if __name__ == "__main__":
    print "start the pipeline"
    
    # TODO:   read parameters (fasta file or (seq,str1,str2))
    #        read ui parameters: treekin-start structure (mfe or oc) or distribution.
    #        max iterations for 2D-Projections.
    
    # input
    seq = "GGGAAUUAUUGUUCCCUGAGAGCGGUAGUUCUC"
    (mfe_struct, mfe) = RNA.fold(seq)
    
    # sampling iterations
    iterations = 1
    
    # input treekin
    startStructureID = 1  # TODO: find the "right" start structure.
    startTime = "0.001"
    endTime = "1e20"
    
    
    # output
    ratesFilePath = "./ratesFile.txt"
    structuresFilePath = "./structuresFile.txt"
    kineticFile = "./treekin-kinetic.txt"
    
    # call samplegenerator which calls RNAxplorer (in matrix2d objects) and executes a clustering (2D or DIANA, etc.) only to select 
    # new structure representatives for samplegeneration.
    structures = samplegenerator.mainloop(seq, mfe_struct, mfe, iterations)

    # call final clusteralg.
    clusters = DIANA.doClustering(structures, 8)
    representatives = selectRepresentatives(clusters)
    
    # prepare representatives and sort them, to make the comparison of the output easier.
    repsAndEnergies = [ (x, RNA.energy_of_struct(seq, x)) for x in representatives]
    repsAndEnergies.sort(key=operator.itemgetter(1, 0))
    
    # rate computation.
    rateMatrix = Rates.arrheniusRates(seq, repsAndEnergies)
    
    # write files for treekin (matrix and structures)
    fw = filewriter.FileWriter()
    fw.writeMatrixToFile(rateMatrix, ratesFilePath)
    fw.writeRepresentativesToFile(seq, repsAndEnergies, structuresFilePath)
    
    # run treekin
    to_call = "cat " + structuresFilePath + " | treekin --p0 " + str(startStructureID) + "=1 --ratesfile " + ratesFilePath + " --t0 " + \
            startTime + " --t8 " + endTime + " -m I > " + kineticFile
    os.system(to_call)
    
    # visualize via webinterface
    

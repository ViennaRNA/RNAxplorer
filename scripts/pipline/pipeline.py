#!/usr/bin/env python

import RNA
import samplegenerator
from clusteralgorithms.diana import DIANA
from rates import Rates
import filewriter
import os
import operator
import sys
import math
from parameters import Parameters, TwoDSamplingParameters

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
    representatives = []
    for c in clusters:
        # select mfe or TODO: centroid structure.
        s = mfeStructure(seq, c)
        representatives.append(s)
    return representatives
    

if __name__ == "__main__":
    print "start the pipeline"
    
    # TODO:   read parameters (fasta file or (seq,str1,str2))
    #        read ui parameters: treekin-start structure (mfe or oc) or distribution.
    #        make command line arguments for all parameters.
    
    # input
    seq = "GGGAAUUAUUGUUCCCUGAGAGCGGUAGUUCUC"
    (mfe_struct, mfe) = RNA.fold(seq)

    # Initial step:
    # start with mfe_struct and open chain as references
    # to retrieve first sample set
    ref_structs = []
    openchain = "." * len(seq)
    ref_structs.append(mfe_struct)
    ref_structs.append(openchain)


    param = Parameters()
    twoDSamplingParameters = TwoDSamplingParameters()
    twoDSamplingParameters.Sequence = seq
    #2D sampling
    twoDSamplingParameters.Reference_one = openchain
    twoDSamplingParameters.Reference_two = mfe_struct
    twoDSamplingParameters.MaxXplorerSamples = 100
    twoDSamplingParameters.SamplingIterations = 10
    param.TwoDSamplingParam = twoDSamplingParameters
    # input treekin
    param.StartTime = 0.001
    param.EndTime = 1e20
    # output
    param.RatesFilePath = "./ratesFile.txt"
    param.StructuresFilePath = "./structuresFile.txt"
    param.KineticFile = "./treekin-kinetic.txt"
    #Diana Threshold (max basepair difference in cluster)
    maxDiameterThreshold = 8
    maxAverageDiameterThreshold = 8
    
    # call samplegenerator which calls RNAxplorer (in matrix2d objects) and executes a clustering (2D or DIANA, etc.) only to select 
    # new structure representatives for samplegeneration.
    structures = samplegenerator.mainloop(twoDSamplingParameters)

    # call final clusteralg.
    clusters = DIANA.doClustering(structures, maxDiameterThreshold, maxAverageDiameterThreshold)
    DIANA.printClusters(clusters)
    representatives = selectRepresentatives(seq, clusters)
    
    # prepare representatives and sort them, to make the comparison of the output easier.
    repsAndEnergies = [ (x, RNA.energy_of_struct(seq, x)) for x in representatives]
    repsAndEnergies.sort(key=operator.itemgetter(1, 0))
    
    # rate computation.
    rateMatrix = Rates.arrheniusRates(seq, repsAndEnergies)
    
    # write files for treekin (matrix and structures)
    fw = filewriter.FileWriter()
    fw.writeMatrixToFile(rateMatrix, param.RatesFilePath)
    fw.writeRepresentativesToFile(seq, repsAndEnergies, param.StructuresFilePath)
    
    # run treekin
    to_call = "cat " + param.StructuresFilePath + " | treekin --p0 " + str(param.StartStructureID) + "=1 --ratesfile " + param.RatesFilePath + " --t0 " + \
            str(param.StartTime) + " --t8 " + str(param.EndTime) + " -m I > " + param.KineticFile
    os.system(to_call)
    
    
    
    

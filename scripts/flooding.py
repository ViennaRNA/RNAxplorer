"""
Script for 2D flooding algorithms and helper functions.
"""

import sys
import operator #for sort with 2 criteria.

def moveSet(klIndex):
    k=klIndex[0]
    l=klIndex[1]
    #Follow the diagonals: (neighbors with bpDist <=2).
    #(Since the k,l-neighborhood is a checkerbord like representation, we have no horizontal or vertical neighbors.)
    return [[k-1,l-1],[k+1,l+1],[k-1,l+1],[k+1,l-1]]

#dictionary for searching the kl-neighbors.
Dictionary={}

#dictionary for searching the kl-entry for a structure.
KLentries={}

def initDictionariesForGradientWalk(landscapeData):
    for entry in landscapeData:
        k=entry[0]
        l=entry[1]
    Dictionary[[k,l]]=entry

    for entry in landscapeData:
        k=entry[0]
        l=entry[1]
        structure=entry[3]
        KLentries[structure]=[k,l]

def gradWalk(structure):
    klIndex = KLentries[structure]
    data=Dictionary[klIndex]
    energy=data[2]
    
    neighborIndices=moveSet(klIndex)
    neighbors = []
    for i in range(0,len(neighborIndices)):
        neighbor=Dictionary[neighborIndices[i]]
        if(neighbor):
            neighbors.append(neighbor)

    if(len(neighbors) == 0):
        return klIndex
    #sort by free energy (and lexicographically to break ties '(' < ')' < '.').
    neighbors = sorted(neighbors,key=operator.itemgetter(2,3))

    newMinKL=[neighbors[0][0],neighbors[0][1]]
    minEnergy=neighbors[0][2]
    neighborStructure=neighbors[0][3]
    
    if(minEnergy < energy | (minEnergy == energy & neighborStructure <= data[3] )):
        return gradWalk(neighborStructure)
    
    return klIndex


def iterativeGradientWalks(landscapeData):
    """"landscapeData=[ [2,5,-3,"...)"],
                    [2,3,-3,".).)"],
                    [0,3,-3,".(.)"],
                    [1,4,-4,".).("],
                    [0,5,-1,".).."]
                  ]"""
    if len(landscapeData) <= 0 | len(landscapeData[0]) < 4:
        print "Error: landscapeData has the wrong structure."
        return []
    
    initDictionariesForGradientWalk(landscapeData)

    #return clusterIndices
    knownMinima=[]
    indexList=[]
    for i in range(0,len(landscapeData)):
        k=landscapeData[i][0]
        l=landscapeData[i][1]
        structure=landscapeData[i][3]

        minKL=gradWalk(structure)
        cID=knownMinima.indexOf(minKL)
    if(cID==-1):
        knownMinima.append(minKL)
        cID=len(knownMinima)-1

        indexList.append([k,l,cID])

    return indexList


def flood2min(landscapeData):
    """
    flood only the two deepest minima up to the height of the first saddle.
    """

    if(landscapeData.length <= 0 | landscapeData[0].length < 4):
        print "Error: landscapeData has the wrong structure."
        return null

    #begin with flooding algorithm
    todoQueue=[]
    doneList=[]
    clusterID=0

    #sort by free energy (and lexicographically to break ties '(' < ')' < '.').
    sortedData = sorted(landscapeData,key=operator.itemgetter(2,3))
    globalMinimum = sortedData[0]
    todoQueue.append(globalMinimum)
    indicesAndClusterIDs=[]
    indicesAndClusterIDs.append([globalMinimum[0],globalMinimum[1],1])
    
    secondMinFound=false
    secondMin=null
    saddleEnergyMaximum=sys.float_info.max
    while(todoQueue.length > 0):  
        currentStructure = todoQueue.shift()
        
        if(doneList.indexOf(currentStructure) == -1):
            doneList.append(currentStructure)
        
            if(currentStructure[2] > saddleEnergyMaximum):
                continue
        
            klIndex = [currentStructure[0],currentStructure[1]]
            neighborIndices = moveSet(klIndex)
            neighbors = []
            for i in range(0, len(neighborIndices)):
                neighbor=Dictionary[neighborIndices[i]]
                if(neighbor):
                    neighbors.append(neighbor)
                
            if len(neighbors) != 0:
                #sort by free energy (and lexicographically to break ties '(' < ')' < '.').
                neighbors = sorted(neighbors,key=operator.itemgetter(2,3))
    
                for neighbor in neighbors:
                    neighborEnergy=neighbor[2]
                    if(neighborEnergy < saddleEnergyMaximum):
                        if(todoQueue.indexOf(neighbor) == -1 & doneList.indexOf(neighbor) == -1):
                            todoQueue.append(neighbor)
                            structure=neighbor[3]
                            minKL=gradWalk(structure)
                            minData=Dictionary[minKL]
                            if(minData[3] == globalMinimum[3]):
                                indicesAndClusterIDs.append([neighbor[0],neighbor[1],1])
    
                            else:
                                #second min found.
                                if(secondMinFound == false):
                                    secondMinFound=true
                                    secondMin=minData
                                    saddleEnergyMaximum=Math.max(neighbor[2],currentStructure[2])
                                    #remove minima with saddleEnergy < saddleEnergyMaximum from the ResultList.
                                    entriesToRemove=[]
                                    for i in range(0,len(indicesAndClusterIDs)):
                                        entry = indicesAndClusterIDs[i]
                                        if(Dictionary[[entry[0],entry[1]]][2] > saddleEnergyMaximum):
                                            entriesToRemove.append(i)
                                        for index in entriesToRemove:
                                            indicesAndClusterIDs.splice(index,1)

                                if(secondMin[2] == minData[2]):
                                    indicesAndClusterIDs.append([neighbor[0],neighbor[1],2])
    
            else:
                console.log("Error: algorithm has stopped, because the cell has zero neighbors!")
                indicesAndClusterIDs.append([currentStructure[0],currentStructure[1],1])
                break
    return indicesAndClusterIDs

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
    return [(k-1,l-1),(k+1,l+1),(k-1,l+1),(k+1,l-1)]

#dictionary for searching the kl-neighbors.
Dictionary={}

#dictionary for searching the kl-entry for a structure.
KLentries={}

def initDictionariesForGradientWalk(landscapeData):
    for entry in landscapeData:
        k=entry[0]
        l=entry[1]
        Dictionary[(k,l)] = entry

    for entry in landscapeData:
        k = entry[0]
        l = entry[1]
        structure = entry[3]
        KLentries[structure] = (k,l)

def gradWalk(structure):
    klIndex = KLentries[structure]

    data = Dictionary[klIndex]
    energy = data[2]
    
    neighborIndices = moveSet(klIndex)
    neighbors = []
    for i in range(0,len(neighborIndices)):
        try:
            neighbor = Dictionary[neighborIndices[i]]
        except KeyError:
            neighbor = False
        if(neighbor):
            neighbors.append(neighbor)

    if(len(neighbors) == 0):
        return klIndex
    #sort by free energy (and lexicographically to break ties '(' < ')' < '.').
    neighbors = sorted(neighbors,key=operator.itemgetter(2,3))

    newMinKL = (neighbors[0][0],neighbors[0][1])
    minEnergy = neighbors[0][2]
    neighborStructure=neighbors[0][3]
    isSmaller = (minEnergy < energy) | ( (minEnergy == energy) & (neighborStructure <= data[3]) )
    if(isSmaller):
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
    knownMinima = []
    indexList = []
    for i in range(0,len(landscapeData)):
        k = landscapeData[i][0]
        l = landscapeData[i][1]
        structure = landscapeData[i][3]

        minKL = gradWalk(structure)
        try:
            cID = knownMinima.index(minKL)
        except ValueError:
            knownMinima.append(minKL)
            cID = len(knownMinima)-1
            
        indexList.append([k,l,cID])

    return indexList


def flood2min(landscapeData):
    """
    flood only the two deepest minima up to the height of the first saddle.
    """

    if( (len(landscapeData) <= 0) | (len(landscapeData[0]) < 4) ):
        print "Error: landscapeData has the wrong structure."
        return []
    
    initDictionariesForGradientWalk(landscapeData)

    #begin with flooding algorithm
    todoQueue = []
    doneList = []
    clusterID = 0

    #sort by free energy (and lexicographically to break ties '(' < ')' < '.').
    sortedData = sorted(landscapeData,key=operator.itemgetter(2,3))
    globalMinimum = sortedData[0]
    todoQueue.append(globalMinimum)
    indicesAndClusterIDs = []
    indicesAndClusterIDs.append([globalMinimum[0],globalMinimum[1],1])
    
    secondMinFound = False
    secondMin = []
    saddleEnergyMaximum=sys.float_info.max
    while(len(todoQueue) > 0):  
        currentStructure = todoQueue.pop()
        try:
            doneList.index(currentStructure)
        except ValueError:
            doneList.append(currentStructure)
        
        if(currentStructure[2] > saddleEnergyMaximum):
            continue
    
        klIndex = (currentStructure[0], currentStructure[1])
        neighborIndices = moveSet(klIndex)
        neighbors = []
        neighbor = []
        for i in neighborIndices:
            try:
                neighbor = Dictionary[i]
            except KeyError:
                neighbor = False
            if(neighbor):
                neighbors.append(neighbor)
            
        if len(neighbors) != 0:
            #sort by free energy (and lexicographically to break ties '(' < ')' < '.').
            neighbors = sorted(neighbors, key=operator.itemgetter(2,3))

            for neighbor in neighbors:
                neighborEnergy = neighbor[2]
                if(neighborEnergy < saddleEnergyMaximum):
                    notInTodo = False
                    notInDone = False
                    try:
                        todoQueue.index(neighbor)
                    except ValueError:
                        notInTodo = True
                    try:
                        doneList.index(neighbor)
                    except ValueError:
                        notInDone = True
                        
                    if(notInDone & notInTodo):
                        todoQueue.append(neighbor)
                        structure = neighbor[3]
                        minKL = gradWalk(structure)
                        minData = Dictionary[minKL]
                        if(minData[3] == globalMinimum[3]):
                            indicesAndClusterIDs.append( [neighbor[0], neighbor[1], 1] )

                        else:
                            #second min found.
                            if(secondMinFound == False):
                                secondMinFound = True
                                secondMin = minData
                                saddleEnergyMaximum = max(neighbor[2],currentStructure[2])
                                #remove minima with saddleEnergy < saddleEnergyMaximum from the ResultList.
                                entriesToRemove = []
                                for i in range(0, len(indicesAndClusterIDs)):
                                    entry = indicesAndClusterIDs[i]
                                    if(Dictionary[(entry[0], entry[1])][2] > saddleEnergyMaximum):
                                        entriesToRemove.append(i)
                                for index in entriesToRemove:
                                    indicesAndClusterIDs.pop(index)

                            if(secondMin[2] == minData[2]):
                                indicesAndClusterIDs.append( [neighbor[0], neighbor[1], 2] )

        else:
            print("Error: algorithm has stopped, because the cell has zero neighbors!")
            #indicesAndClusterIDs.append( [currentStructure[0], currentStructure[1], 1] )
            break

    return indicesAndClusterIDs

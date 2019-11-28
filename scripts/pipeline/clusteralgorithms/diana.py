#!/usr/bin/python3
"""
Diana Clustering module for clustering RNA structures according to their basepairdistance.
The cluster distance function computes the distance between the centroids of the clusters.

! This algorithm may present different results than the algorithm in R. 
The reason is that the maximal diameter (criterion for selecting cluster to split) is not unique !
A second reason is that the object with maximal average distance is not unique. 

"""

import sys, math, RNA, numpy, argparse, re

class Cluster:
    """
    diana internal tree-like cluster structure.
    """
    def __init__(self):
        self.structures = []
        self.representative = ""
        self.childNodes = []

def computeBasePairDistanceMatrix(structs):
    matrix = numpy.zeros((len(structs), len(structs)), dtype=numpy.int)
    for i in range(len(structs)):
        for j in range(i + 1, len(structs)):
            dist = RNA.bp_distance(structs[i], structs[j])
            matrix[i][j] = dist
            matrix[j][i] = dist
    return matrix

def averageDissimilarity(index, cluster, dissMatrix):
    """
    computes the average dissimilarity
    index is a structure index, which is conform the the index in the dissMatrix.
    cluster is a list of structure indices.
    """
    avgDiss = 0.0
    subtract = 0.0
    if index in cluster:
        subtract = 1.0
    for i in cluster:
        avgDiss += float(dissMatrix[index][i])
    NumberOfElementsWithoutIndexElement = float(len(cluster) - subtract)
    if NumberOfElementsWithoutIndexElement > 0.0:
        avgDiss = avgDiss / NumberOfElementsWithoutIndexElement
    else:
        avgDiss = 0.0
    return avgDiss
    
def objectIndexWithHighestDissimilarity(cluster, dissMatrix):
    maxAvgDiss = 0.0
    maxIndex = -1
    for i1 in cluster:
        avgDiss = averageDissimilarity(i1, cluster, dissMatrix)
        if (avgDiss > maxAvgDiss):
            maxAvgDiss = avgDiss
            maxIndex = i1
    return maxIndex     
    
def diameter(cluster, dissMatrix):
    """
    computes the largest basepairdistance between all pairs of structures in the cluster.
    """
    maxDist = 0
    for i in range(len(cluster)):
        for j in range(i + 1, len(cluster)):
            dist = dissMatrix[cluster[i]][cluster[j]]
            if dist > maxDist:
                maxDist = dist
    return maxDist

def WmAll(c, dissMatrix):
    """Computes the 'within-ness' measure for a cluster, with all pairwise comparisons.

    Args:
        c (list): List of secondary structures.

    Returns:
        set: The within-ness measure.
    """
    distanceSum = 0.0
    for i in range(0, len(c)):
        for j in range(i + 1, len(c)):
            distanceSum += math.pow(dissMatrix[c[i]][c[j]], 2)
    return distanceSum

def sum_distance_in_cluster(c, dissMatrix):
    sum_dist = 0.0
    for i in range(0, len(c)):
        for j in range(i + 1, len(c)):
            sum_dist += abs(dissMatrix[c[i]][c[j]])
    return sum_dist

def mean_distance_in_cluster(c, dissMatrix):
    mean_dist = sum_distance_in_cluster(c, dissMatrix)
    elements = ((len(c) * (len(c) -1)) / 2.0)
    if elements <= 0:
        return mean_dist
    else:
        mean_dist = mean_dist / elements
    return mean_dist

def mean_distance_in_matrix(dissMatrix, matrix_size):
    avg_dist = 0
    for i in range(matrix_size):
        for j in range(i+1, matrix_size):
            avg_dist += abs(dissMatrix[i][j])
    n_dist = (float(matrix_size) * float(matrix_size-1))/2.0
    avg_dist = avg_dist / n_dist
    return avg_dist
    

def varianceRatio(parentCluster, clusterA, clusterB, dissMatrix):
    """
    compute the variance of the old clusters and the new cluster.
    Args:
        cl1 (list): Former clustering, list of lists of secondary structures.
        cl2 (list): Current clustering, list of lists of secondary structures.
    """
    
    wmA = WmAll(clusterA.structures, dissMatrix)
    wmB = WmAll(clusterB.structures, dissMatrix)
    wmParent = WmAll(parentCluster.structures, dissMatrix)
    
    withinRatio = (wmA + wmB) / wmParent
    randomRatio = (pow(len(clusterA.structures), 2.) + pow(len(clusterB.structures), 2.)) / (pow(len(parentCluster.structures), 2.))
    
    return withinRatio, randomRatio

def pooled_within_group_sum_of_squares(clusters, dissMatrix):
    pwgss = 0
    if clusters == None:
        return pwgss
    for c in clusters:
        pwgss += sum_distance_in_cluster(c, dissMatrix) / float(len(c))
        """
        faster alternatice is could be sum_i_N(pow(bpdist(centroid_structure(c), structure(c_i)),2))
        """
    return pwgss

def pooled_between_group_sum_of_squares(clusters, dissMatrix, matrix_size):
    pbgss = 0
    if clusters == None:
        return pbgss
    avg_dist = mean_distance_in_matrix(dissMatrix, matrix_size)
            
    for c in clusters:
        avg_dist_c = mean_distance_in_cluster(c, dissMatrix)
        pbgss +=  abs(avg_dist_c - avg_dist) * float(len(c))
        """
        faster alternatice is could be sum_i_N(pow(bpdist(centroid_structure(dataset), centroid_structure(c_i)),2)) /float(len(c))
        """
    return pbgss

def calinski_harabasz_index(clusters, dissMatrix, matrix_size):
    ch = 0
    pbgss = pooled_between_group_sum_of_squares(clusters, dissMatrix, matrix_size)
    pwgss = pooled_within_group_sum_of_squares(clusters, dissMatrix)
    number_of_clusters = len(clusters)
    number_of_structures = matrix_size
    ch =  ((pbgss) / (pwgss)) * ((number_of_structures - number_of_clusters) / float(number_of_clusters -1))
    return ch
    

def maxDiameterInLeafNodes(clusterTree, dissMatrix):
    maxDiameter = -1
    maxCluster = None
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            dia, cluster = maxDiameterInLeafNodes(cn, dissMatrix)
            if dia > maxDiameter:
                maxDiameter = dia
                maxCluster = cluster
    else:
        maxDiameter = diameter(clusterTree.structures, dissMatrix)
        maxCluster = clusterTree
    
    return maxDiameter, maxCluster

def averageDiameterInLeafNodes(clusterTree, dissMatrix):
    """
    also known as Divisive Coefficient (DC).
    """
    sumDiameters = 0
    numberOfClusters = 0
    stackChildNodes = []
    stackChildNodes.append(clusterTree)
    while(len(stackChildNodes) > 0):
        cn = stackChildNodes.pop()
        if len(cn.childNodes) == 0:
            #is leaf node/cluster --> sum
            sumDiameters += diameter(cn.structures, dissMatrix)
            numberOfClusters += 1
        else:    
            for c in cn.childNodes:
                stackChildNodes.append(c)
                
    averageDiameter = sumDiameters / float(numberOfClusters)
    return averageDiameter

def createClusterTree(c_root, dissMatrix, maxDiameterThreshold, maxAverageDiameterThreshold, structs, do_ch_first_local_min = False, hierarchy=1):
    """
    The core of the diana algorithm (recursive function).
    c_root = the rootnode of the clusterTree. It contains the main cluster as childnode.
    """
    # 1. select cluster with the largest diameter from all leafnodes.
    dc = averageDiameterInLeafNodes(c_root, dissMatrix)
    if dc < maxAverageDiameterThreshold:
        return
    if maxDiameterThreshold >= 0:
        maxDiameter, c_m = maxDiameterInLeafNodes(c_root, dissMatrix)
        if maxDiameter <= maxDiameterThreshold:
            return      
            
    if c_m == None:
        return
    
    if len(c_m.structures) > 1:
        hierarchy+=1
        # 2. object with highest dissimilarity to all others defines the new cluster (c_newA).
        o = objectIndexWithHighestDissimilarity(c_m.structures, dissMatrix)
        c_newA = Cluster()
        c_newA.structures.append(o)
        # B contains all other elements
        c_newB = Cluster()
        c_newB.structures.extend([x for x in c_m.structures if x != o])
        
        # 3. select best elements for each cluster (move closest elements from B to A).
        # for each object outside the splinter group:
        positiveDiExists = True
        while(positiveDiExists):
            positiveDiExists = False
            largestD_i = 0.0
            bestElement = None
            for j in c_newB.structures:
                d_i = averageDissimilarity(j, c_newB.structures, dissMatrix) - averageDissimilarity(j, c_newA.structures, dissMatrix)
                if (d_i > largestD_i):
                    largestD_i = d_i
                    bestElement = j
            if (bestElement != None) & (len(c_newB.structures) > 1):
                positiveDiExists = True
                c_newA.structures.append(bestElement)
                c_newB.structures.remove(bestElement)
        
        # ratios = varianceRatio(c_m, c_newA, c_newB, dissMatrix)
        # criterion = ratios[0] < ratios[1]
        # print ratios[0], ratios[1]
        # if(criterion):
        c_m.childNodes.append(c_newB)
        c_m.childNodes.append(c_newA)

        if do_ch_first_local_min != False:
            if do_ch_first_local_min == True:
                do_ch_first_local_min = sys.float_info.max
            cis = convertTreeToListOfClusters_indices(c_root, structs, [])
            #print(cis, len(cis))
            if len(cis) > 1:
                ch = calinski_harabasz_index(cis, dissMatrix, len(structs))
                if ch > do_ch_first_local_min:
                    if len(c_m.childNodes) > 0:
                        c_m.childNodes = []
                        print(hierarchy, "ch", "{:10.3f}".format(ch))
                    return # ch index > last ch index --> first local min --> break
                print(hierarchy, "ch", "{:10.3f}".format(ch))
                do_ch_first_local_min = ch

        # build the subtrees for each child
        createClusterTree(c_root, dissMatrix, maxDiameterThreshold, maxAverageDiameterThreshold, structs, do_ch_first_local_min, hierarchy)


def convertTreeToListOfClusters(clusterTree, structs, clusterList=[]):
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            convertTreeToListOfClusters(cn, structs, clusterList)
    else:
        cluster = [ structs[c] for c in clusterTree.structures]
        clusterList.append(cluster)
    return clusterList

def convertTreeToListOfClusters_indices(clusterTree, structs, clusterList=[]):
    if len(clusterTree.childNodes) > 0:
        for cn in clusterTree.childNodes:
            convertTreeToListOfClusters_indices(cn, structs, clusterList)
    else:
        cluster = [ c for c in clusterTree.structures]
        clusterList.append(cluster)
    return clusterList 
         
class DIANA:
    @staticmethod
    def doClustering(structs, maxDiameter, maxAverageDiameter, do_ch_first_local_min = False):
        """
        Computes the DIANA clustering
        
        Args:
            maxDiameter = threshold for clustering
            maxAverageDiameter = "
        """
    
        structs = list(structs)
        # start DIANA
        # 1. initialization
        dissMatrix = computeBasePairDistanceMatrix(structs)
        c_root = Cluster()
        c_root.structures.extend([x for x in range(0, len(structs))])
        
        # start the real DIANA algorithm.
        createClusterTree(c_root, dissMatrix, maxDiameter, maxAverageDiameter, structs, do_ch_first_local_min)
        # end DIANA
              
        clusters = convertTreeToListOfClusters(c_root, structs, [])
        
        #cis = convertTreeToListOfClusters_indices(c_root, structs)
        #if len(cis) > 1:
        #    ch = calinski_harabasz_index(cis, dissMatrix, len(structs))
        #    print("ch", ch)
        #else:
        #    print("ch", 0)
        return clusters
    
    @staticmethod
    def printClusters(clusters):
        """
        print a list of lists.
        """
        if clusters == None:
            return
        cid = 0
        for c in clusters:
            cid +=1
            print("ClusterID:",cid)
            for s in c:
                print(s)

def parseFile(fpath):
    """
    Parses a FASTA-ish file and extract a list of secondary structures.

    Args:
        fpath (string): Path to file in FASTA-ish notation, featuring a list
        of Vienna-formatted secondary structures.

    Returns:
        list: List of secondary structures, each represented as a pair `(v,bps)`,
        where `v` is the Vienna notation and `bps` is a set of base-pairs.
    """
    res = []
    for l in open(fpath):
        if not l.startswith("#") and not l.startswith(";") and not l.startswith(">"):
            match = re.search('([\.\(\)]+)', l)
            if match:
                sequence = match.group(1)
                res.append(sequence.strip())
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Implememntation of the DIANA divisive clustering for RNA secondary structures.')
    parser.add_argument("-f", "--file", type=str, required=True, help="Fasta file with secondary structures")
    parser.add_argument("-d", "--diameter-threshold", type=int, default=0, required=False, help="Cluster diameter threshold (max bp distance within a cluster)")
    parser.add_argument("-m", "--average-diameter-threshold", type=float, default=0, required=False, help="Average diameter threshold")
    parser.add_argument("-c", "--calinski-harabasz-threshold", action='store_true', default=False, required=False, help="Abort if the Calinski-Harabasz-Index reaches the first local minimum")
    args = parser.parse_args()
    structureFileName = args.file #sys.argv[1]
    structs = parseFile(structureFileName)
    maxDiameter = args.diameter_threshold
    maxAverageDiameter = args.average_diameter_threshold
    d = DIANA()
    clusters = d.doClustering(structs, maxDiameter, maxAverageDiameter, args.calinski_harabasz_threshold)
    d.printClusters(clusters)
    









from tinyRNAClustering import *
import sys,os

RNAFOLD_SAMPLING = "rnasubopt -p 10000 -s <"
R_EXEC = "RScript"

def unique(structs):
    res = {}
    for s,ss in structs:
        res[s] = ss
    return [(x,res[x]) for x in res]

def generateStructs(seq,outpath):
    outfile = open("tmp.faa","w")
    outfile.write(seq+"\n")
    outfile.close()
    os.system(RNAFOLD_SAMPLING+ " tmp.faa > "+outpath)
    return unique(parseFile(outpath))
    

def computeDisMatrix(structs):
    return [[computeDistance(ss1,ss2) for j,(s2,ss2) in enumerate(structs)] for i,(s1,ss1) in enumerate(structs)]
    
def exportDisMatrix(m,fp):
    f = open(fp,"w")
    for i in range(len(m)):
        for j in range(len(m[i])):
            if j!=0:
                f.write(",") 
            f.write("%s"%m[i][j])            
        f.write("\n")
    f.close()

def runDIANA(matPath):
    os.system(R_EXEC+ " ..\R\pw_distance.R "+matPath)
    return matPath+".clusters.csv"


def createClusters(seq):
    structsPath = "structs.fa"
    disMatPath = "distMat.dat"
    structs = generateStructs(seq,structsPath)
    mat = computeDisMatrix(structs)
    exportDisMatrix(mat,disMatPath)
    clustersPath = runDIANA(disMatPath)
    doClustering(structs,clustersPath)
    

if __name__=="__main__":
    createClusters(sys.argv[1])

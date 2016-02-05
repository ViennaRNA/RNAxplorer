"""
Script for 2D flooding algorithms and helper functions.
"""

import operator  # for sort with 2 criteria.

"""
Objects of the RNADfoldColumn class represent a column in a
2D projection with respect to some reference structures

In essence, this just collects a number of secondary structures
in dot-bracket notation, and stores its mfe representative. Thus,
it can be used for any projection, no matter the number of dimensions
"""
class RNA2DfoldColumn(object):

    def __init__(self, structure, e, idx):
        self.structures = []
        self.mfe = 10000.0
        self.mfe_struct = ""

        self.addStructure(structure, e, idx)

    def addStructure(self, s, e, idx):
        self.structures.append((s, idx))
        if self.mfe > e:
            self.mfe = e
            self.mfe_struct = s

class WatershedFlooder:
    @staticmethod
    def doflooding(landscapeData):
        """"landscapeData=[ [2,5,-3,"...)"],
                        [2,3,-3,".).)"],
                        [0,3,-3,".(.)"],
                        [1,4,-4,".).("],
                        [0,5,-1,".).."]
                      ]"""
        if len(landscapeData) <= 0 | len(landscapeData[0]) < 4:
            print "Error: landscapeData has the wrong structure."
            return []
    
        # put data into a 2D grid (actually a dict with (k,l) keys)
        grid2D = {}
        for i in range(0, len(landscapeData)):
            k = landscapeData[i][0]
            l = landscapeData[i][1]
            e = landscapeData[i][2]
            s = landscapeData[i][3]
            if (k, l) in grid2D:
                grid2D[(k, l)].addStructure(s, e, i)
            else:
                grid2D[(k, l)] = RNA2DfoldColumn(s, e, i)
    
        # create a to-do list of columns to process
        todo = []
        for k in grid2D.keys():
            e = grid2D[k].mfe
            todo.append((k, e))
    
        # sort the to-do list by free energy
        todo.sort(key=operator.itemgetter(1))
    
        # list of already processed columns
        done = []
    
        basins = []
    
        for i in range(0, len(todo)):
            basin = []
            todo_stack = []
    
            ((k, l), e) = todo[i]
            
            todo_stack.append((k, l))
    
            while len(todo_stack) > 0:
                # cell is a (k,l) tuple
                cell = todo_stack.pop()
    
                if not cell in done:
                    basin.append(cell)
    
                    cell_e = grid2D[cell].mfe
    
                    # determine neighbor grid cells of current cell
                    neighbors = []
                    nn = [(k - 1, l - 1), (k + 1, l - 1), (k - 1, l + 1), (k + 1, l + 1)]
                    for n in nn:
                        if n in grid2D:
                            neighbors.append(n)
    
                    for n in neighbors:
                        n_e = grid2D[n].mfe
                        if n_e >= cell_e and n not in basin:
                            todo_stack.append(n)
    
            if len(basin) > 0:
                basins.append(basin)
    
            for b in basin:
                done.append(b)
    
        # now collect all the structures again
        bbb = []
    
        def combineBasins(b):
            c = []
            for n in b:
                map(lambda v: c.append(v), grid2D[n].structures)
            bbb.append(c)
    
        map(combineBasins, basins)
    
    
        return bbb

#!/usr/bin/python3

import argparse
import re
import RNA
import operator
import sys
from Bio import Phylo
import pylab
from math import *
from matplotlib.pyplot import figure
from matplotlib import pyplot
from multiprocess.pool import Pool
import numpy as np


class StructureEnergy():
    Index = None
    Structure = None
    Energy = None
    def __init__(self, index, structure, energy):
        self.Index = index
        self.Structure = structure
        self.Energy = energy


"""
- assumes sorted indices in barriers file with one based index!
- structure_list has one based indices
"""   
def read_structure(file):
    sequence = None
    structure_list = []
    min_energy = sys.maxsize
    with open(file, 'r') as f:
        new_index = 0
        for line in f:
            match = re.match("\s*(\d+)?\s*([\.\(\)]+)\s*(\-?\d+\.?\d*)", line)
            if match:
                index = match.group(1)
                structure = match.group(2)
                energy = match.group(3)
                energy = float(energy)
                if energy < min_energy:
                    min_energy = energy
                try:
                    index = int(index)
                except Exception as e:
                    new_index += 1
                    index = new_index
                structure_object = StructureEnergy(index, structure, energy)
                structure_list.append(structure_object)
                #print(index, structure, energy)
            else:
                match = re.match("\s*([ACGUTNIP]+)", line)
                if match:
                    sequence = match.group(1)
    return sequence, structure_list, min_energy

def read_pure_structure_list(file):
    structure_list = []
    with open(file, 'r') as f:
        new_index = 0
        for line in f:
            match = re.match("([\.\(\)]+)", line)
            if match:
                structure = match.group(1)
                energy = None
                new_index += 1
                index = new_index
                structure_object = StructureEnergy(index, structure, energy)
                structure_list.append(structure_object)
    return structure_list


"""
- assumes sorted indices in barriers file with one based index!
- structure_list has one based indices
- saddle_list has zero based indices!
"""
def read_bar_saddle_file(file, filter):
    structure_count = 0
    sequence = None
    structure_list = []
    saddle_list = []
    min_energy = sys.maxsize
    with open(file, 'r') as f:
        new_index = 0
        for line in f:
            match = re.match("\s*(\d+)\s*([\.\(\)]+)\s*(\-?\d+\.?\d*)\s*(\d+)\s*(\-?\d+\.?\d*)\s*([\.\(\)\~]+)?", line)
            if match:
                structure_count +=1
                if filter.MaxMinima != None and structure_count > filter.MaxMinima:
                    break
                
                index = match.group(1)
                structure = match.group(2)
                energy = match.group(3)
                energy = float(energy)
                structure_basin_index = match.group(4)
                structure_basin_index = int(structure_basin_index)
                saddle_energy_distance = match.group(5)
                saddle_energy_distance = float(saddle_energy_distance)
                saddle_energy_kcal = saddle_energy_distance + energy
                saddle_structure = match.group(6)
                #print(index, structure, energy, structure_basin_index,saddle_energy_distance, saddle_structure)
                if energy < min_energy:
                    min_energy = energy
                try:
                    index = int(index)
                except Exception as e:
                    new_index += 1
                    index = new_index
                structure_object = StructureEnergy(index, structure, energy)
                structure_list.append(structure_object)
                #print(index, structure, energy, saddle_energy_kcal)
                saddle_list.append((index-1, structure_basin_index-1, saddle_energy_kcal))
            else:
                match = re.match("\s*([ACGUTNIP]+)", line)
                if match:
                    sequence = match.group(1)
    #print(saddle_list)
    return sequence, structure_list, min_energy, saddle_list

def find_saddle(s1, s2, index_1, index_2):
    fold_compound = RNA.fold_compound(sequence)
    saddle_energy_dcal = fold_compound.path_findpath_saddle(s1, s2)
    saddle_energy_kcal = saddle_energy_dcal / 100.0
    return (index_1, index_2, saddle_energy_kcal)

def connect_structures_find_saddles(sequence, structure_list):
    pairs = {}
    #fc = RNA.fold_compound(sequence)
    fp_pool = Pool(Max_Threads)
    res_list=[]
    for i, se_1 in enumerate(structure_list):
        for j in range(i+1, len(structure_list)):
            se_2 = structure_list[j]
            a = fp_pool.apply_async(find_saddle, args=(se_1.Structure, se_2.Structure, i, j))
            res_list.append(a)
            #saddle_energy_dcal = fc.path_findpath_saddle(se_1.Structure, se_2.Structure)
            #saddle_energy_kcal = saddle_energy_dcal / 100.0
            #pairs[(i,j)] = saddle_energy_kcal
            #pairs[(j,i)] = saddle_energy_kcal
    fp_pool.close()
    
    for a in res_list:
        i,j, saddle_energy_kcal = a.get()
        pairs[(i,j)] = saddle_energy_kcal
        pairs[(j,i)] = saddle_energy_kcal
        
        
    # get lowest saddle for each structure that ends in a structure with lower energy than the first structure.
    minimal_saddle_list = []
    for i in range(0, len(structure_list)):
        se_1 = structure_list[i]
        min_saddle_energy = sys.maxsize
        tree_neighbor = None
        for j in range(0, len(structure_list)):
            if i == j:
                continue
            se_2 = structure_list[j]
            saddle_energy = pairs[(i,j)]
            if saddle_energy <= min_saddle_energy and se_2.Energy < se_1.Energy:
                min_saddle_energy = saddle_energy
                tree_neighbor = j
        if tree_neighbor == None: # it could be the root.
            tree_neighbor = -1
        minimal_saddle_list.append((i, tree_neighbor, min_saddle_energy))
    return minimal_saddle_list


class Node():
    Children = []
    Reachable_Leaf_Indices = None
    Parent = None
    Branch_Length = None
    Saddle_Energy = None
    Energy = None
    Index = None
    def __init__(self, index, energy, saddle_energy, parent = None, children = None):
        self.Index = index
        self.Energy = energy
        self.Saddle_Energy = saddle_energy
        self.Branch_Length = 0.0
        if saddle_energy != None and energy != None:
            self.Branch_Length = self.Saddle_Energy - self.Energy
        self.Parent = parent
        self.Children = []
        if type(children) == type([]):
            self.Children = children
            
        self.Reachable_Leaf_Indices = set()
        if index != None:
            self.Reachable_Leaf_Indices |= set([index])
    
    def extend_branch_length(self,additional_length):
        self.Branch_Length += additional_length
    
    """
    return true if the subtree defined by this node contains the index
    """
    def contains_node_index(self, index):
        has_index = False
        #print(self.Index, index, "comp")
        if self.Index == index:
            has_index = True
        else:
            if len(self.Children) > 0:
                for c in self.Children:
                    #print("ch", c.Index, index)
                    if c.contains_node_index(index):
                        has_index = True
                        break
        return has_index
    
    
    def find_node_index(self, index):
        result_node = None
        if self.Index == index:
            result_node = self
        else:
            if len(self.Children) > 0:
                for c in self.Children:
                    res =  c.find_node_index(index)
                    if res != None:
                        result_node = res
                        break
        return result_node

    def get_leafs(self):
        leafs = set()
        if self.Index != None and len(self.Children) == 0:
            leafs.add(self.Index)
        else:
            if len(self.Children) > 0:
                for c in self.Children:
                    leafs |= c.get_leafs()
        return leafs
    

def create_barrier_tree(minimal_saddle_list, structure_list):
    n_structures = len(structure_list)
    saddle_matrix = []

    saddle_matrix = np.empty((n_structures,n_structures,)) #sparse.lil_matrix((n_structures, n_structures))
    saddle_matrix.fill(np.nan)
    
    clusters = []
    has_updated_length_for_ids = []
    clustered_ids = set()
    previous_saddle = None

    for s in minimal_saddle_list:
        id_from = s[0]
        id_to = s[1]

        if id_to == -1: # ignore unconnected root edge
            continue
        
        saddle_energy = s[2]
        
        saddle_difference = 0
        if previous_saddle != None:
            saddle_difference = saddle_energy - previous_saddle
        previous_saddle = saddle_energy
        #print(id_from, id_to)
        if id_from in clustered_ids and id_to in clustered_ids:
            #update lengths of two branches, if not already done (use lowest saddle of sorted saddle list)
            id_pair = set([id_from, id_to])
            if not id_pair in has_updated_length_for_ids:
                #get branches where i is in one and j in the other
                c_id_to = None
                c_id_from = None
                
                for c_i, c in enumerate(clusters):
                    if c.contains_node_index(id_from):
                        c_id_from = c_i
                        break
                for c_i, c in enumerate(clusters):
                    if c.contains_node_index(id_to):
                        c_id_to = c_i
                        break
                if c_id_from == None or c_id_to == None:
                    print("Error: handled node is not in tree.")
                    exit()
                if c_id_from != c_id_to:
                    #print("update lengths", id_from, id_to)
                    # update lengths only if they are not in the same cluster
                    clusters[c_id_from].extend_branch_length(saddle_difference)
                    clusters[c_id_from].Saddle_Energy = saddle_energy
                    clusters[c_id_to].extend_branch_length(saddle_difference)
                    clusters[c_id_to].Saddle_Energy = saddle_energy
                    # and merge them!
                    
                    parent = Node(None, None, saddle_energy, None, children = [clusters[c_id_from], clusters[c_id_to]])
                    clusters[c_id_from].Parent = parent
                    clusters[c_id_to].Parent = parent
                    
                    set_from = clusters[c_id_from].Reachable_Leaf_Indices
                    set_to = clusters[c_id_to].Reachable_Leaf_Indices 
                    for leaf_index in set_from:
                        for leaf_index_2 in set_to:
                            if np.isnan(saddle_matrix[leaf_index_2, leaf_index]):
                                saddle_matrix[leaf_index_2, leaf_index] = saddle_energy
                                saddle_matrix[leaf_index, leaf_index_2] = saddle_energy
                    
                    
                    parent.Reachable_Leaf_Indices |= clusters[c_id_from].Reachable_Leaf_Indices
                    parent.Reachable_Leaf_Indices |= clusters[c_id_to].Reachable_Leaf_Indices
                    
                    clusters.pop(max(c_id_from, c_id_to))
                    clusters.pop(min(c_id_from, c_id_to))
                    clusters.append(parent)
                else:
                    # update lengths in same cluster
                    clusters[c_id_from].extend_branch_length(saddle_difference)
                    
                has_updated_length_for_ids.append(id_pair)
            continue
        
        if not id_from in clustered_ids and not id_to in clustered_ids:
            #print("new cluster", id_from, id_to)
            # create new cluster with both ids
            energy_from = structure_list[id_from].Energy
            energy_to = structure_list[id_to].Energy
            node_from = Node(id_from, energy_from, saddle_energy, parent = None, children = None)
            node_to = Node(id_to, energy_to, saddle_energy, parent = None, children = None)
            
            index = None
            energy = None
            #saddle_energy = None
            parent = Node(index, energy, saddle_energy, parent = None, children = [node_from, node_to])
            parent.Reachable_Leaf_Indices |= set([id_from, id_to])
            node_from.Parent = parent
            node_to.Parent = parent
            
            for c in clusters:
                #c.extend_branch_length(saddle_difference)
                diff = saddle_energy - c.Saddle_Energy
                c.Branch_Length += diff
                c.Saddle_Energy = saddle_energy
            
            clusters.append(parent)
                
            clustered_ids.add(id_from)
            clustered_ids.add(id_to)
            
            saddle_matrix[id_from, id_to] = saddle_energy
            saddle_matrix[id_to, id_from] = saddle_energy
            continue
        
        if id_from in clustered_ids and not id_to in clustered_ids:
            #print("add cluster", id_to)
            # create cluster for id to and merge it to cluster with id_from
            c_from_id = None
            for c_i, c in enumerate(clusters):
                if c.contains_node_index(id_from):
                    c_from_id = c_i
                    break
            node_from = clusters.pop(c_from_id)
            diff = saddle_energy - node_from.Saddle_Energy
            node_from.Branch_Length += diff
            node_from.Saddle_Energy = saddle_energy

            
            energy_to = structure_list[id_to].Energy
            node_to = Node(id_to, energy_to, saddle_energy, parent = None, children = None)
            
            index = None
            energy = None
            #saddle_energy = None
            parent = Node(index, energy, saddle_energy, parent = None, children = [node_from, node_to])
            parent.Reachable_Leaf_Indices |= node_from.Reachable_Leaf_Indices
            parent.Reachable_Leaf_Indices |= node_to.Reachable_Leaf_Indices
            
            node_from.Parent = parent
            node_to.Parent = parent
            
            for c in clusters:
                #c.extend_branch_length(saddle_difference)
                diff = saddle_energy - c.Saddle_Energy
                c.Branch_Length = +diff
                c.Saddle_Energy = saddle_energy
                
            clusters.append(parent)
            
            clustered_ids.add(id_to)
        
            for leaf_index in node_from.Reachable_Leaf_Indices:
                if np.isnan(saddle_matrix[id_to, leaf_index]):
                    saddle_matrix[id_to, leaf_index] = saddle_energy
                    saddle_matrix[leaf_index, id_to] = saddle_energy
            
            continue
            
        if not id_from in clustered_ids and id_to in clustered_ids:
            #print("add cluster", id_from)
            # create cluster for id_from and merge it to id_to
            c_to_id = None
            for c_i, c in enumerate(clusters):
                if c.contains_node_index(id_to):
                    c_to_id = c_i
                    break
            node_to = clusters.pop(c_to_id)
            diff = saddle_energy - node_to.Saddle_Energy
            node_to.Branch_Length += diff
            node_to.Saddle_Energy = saddle_energy

            energy_from = structure_list[id_from].Energy
            node_from = Node(id_from, energy_from, saddle_energy, parent = None, children = None)
            
            index = None
            energy = None
            #saddle_energy = None
            parent = Node(index, energy, saddle_energy, parent = None, children = [node_from, node_to])
            parent.Reachable_Leaf_Indices |= node_from.Reachable_Leaf_Indices
            parent.Reachable_Leaf_Indices |= node_to.Reachable_Leaf_Indices
            
            node_from.Parent = parent
            node_to.Parent = parent
            
            for c in clusters:
                #c.extend_branch_length(saddle_difference)
                diff = saddle_energy - c.Saddle_Energy
                c.Branch_Length += diff
                c.Saddle_Energy = saddle_energy
            
            clusters.append(parent)
            
            clustered_ids.add(id_from)
            
            for leaf_index in node_to.Reachable_Leaf_Indices:
                if np.isnan(saddle_matrix[id_from, leaf_index]):
                    saddle_matrix[id_from, leaf_index] = saddle_energy
                    saddle_matrix[leaf_index, id_from] = saddle_energy
            continue
    
    return clusters, saddle_matrix


def newick_string_builder(tree):
    if type(tree) != Node:
        return ""
    sub_tree_string ="("
    if tree.Children != None and len(tree.Children) > 0:
        sub_tree_string += "("
        for c in tree.Children:
            sub_tree_string += newick_string_builder(c)
            sub_tree_string += ","
        if sub_tree_string[-1] == ",":
            sub_tree_string = sub_tree_string[:-1]
        sub_tree_string += ")"
    else:
        # it is a leaf node --> add index
        sub_tree_string += str(tree.Index+1)
    sub_tree_string += ":" + "{:.2f}".format(tree.Branch_Length) + ")"
    return  sub_tree_string
    

def filter_tree_min_height(tree, minh):
    """
    TODO: implement this.
    """
    return tree
  
def create_newick_tree_string(minimal_saddle_list, structure_list, filter):
    #apply max Energy filter
    filtered_saddle_list = []
    if filter.MaxEnergy != None:
        for s in minimal_saddle_list:
            if s[2] < filter.MaxEnergy:
                filtered_saddle_list.append(s)
    else:
        filtered_saddle_list = minimal_saddle_list
    
    # create barriers tree (could be one or many if it is unconnected)
    clusters, saddle_matrix = create_barrier_tree(filtered_saddle_list, structure_list)

    filtered_clusters = []
    if filter.MinHeight != None:
        for c in clusters:
            minh_tree = filter_tree_min_height(c, filter.MinHeight)
            filtered_clusters.append(minh_tree)
    else:
        filtered_clusters = clusters
    
    if len(filtered_clusters) <= 0:
        print("Error: no tree was created!")
    if len(filtered_clusters) > 1:
        print("Error: the landscape is not connected! We have a forest!")

    barriers_tree = None
    for i in range(len(filtered_clusters)):
        if filtered_clusters[i].contains_node_index(0):
            # if it contains the mfe basin (with index 0)
            barriers_tree = filtered_clusters.pop(i)
            break

    newick_string = ""
    if barriers_tree != None:
        newick_string = newick_string_builder(barriers_tree) + ";"
        print("tree to plot:")
        print(newick_string)
    else:
        print("Error: the mfe structure is not connected to the tree!")
    
    if len(filtered_clusters) > 0:
        print("other trees")
        for c in filtered_clusters:
            other_tree = newick_string_builder(c) + ";"
            print(other_tree)
            if newick_string == "":
                newick_string = other_tree

    return newick_string, filtered_saddle_list, barriers_tree, saddle_matrix


def print_newick_tree(newick_string, output_file, max_saddle_energy, min_energy, ids_to_color):
    n = Phylo.NewickIO.StringIO(newick_string)
    tree = Phylo.read(n, 'newick')
    """
    mrca = tree.common_ancestor({"name": "1"}, {"name": "2"})
    mrca.color = "salmon"
    #tree.find_clades({'name': '1'})[0].color = 'red'
    leafs = tree.get_terminals()
    for l in leafs:
        l.color = 'red'
    """
    
    # create color hash map for leafs that get a different color
    minima_colors = {}
    for id in ids_to_color:
        minima_colors[str(id)] = 'r'

    # print longest edges at first
    tree.ladderize(reverse=True)

    # inverse the x axis by changing the x-tick-labels (because edges weights in the tree are positive)
    # we could also invert the axis and use negative distances in the tree, but this does not work well with Phylo.draw.
    offset = ceil(max_saddle_energy) - max_saddle_energy
    my_labels = [x for x in range(int(ceil(max_saddle_energy)),int(floor(min_energy))-2, -1)]
    my_ticks = [ x - offset for x in range(len(my_labels))]
    
    pyplot.rcParams["font.size"] = 25
    pyplot.rcParams["lines.linewidth"] = 3
    figure(figsize=(50,30))
    my_axes = pylab.axes(xticks = my_ticks, xticklabels = my_labels)
    Phylo.draw(tree, do_show=False, label_colors=minima_colors, axes=my_axes, branch_labels=lambda c: "{:.2f}".format(float(c.branch_length)) if c.branch_length != None else None)
    #Phylo.draw(tree, do_show=False, label_colors={'1':'r', '2':'g', '3':'g'}, axes=my_axes)
    #pylab.ylim((0,n_minima))
    pylab.xlim((min(my_ticks),max(my_ticks)))
    #my_axes.invert_xaxis()
    pylab.yticks([])
    pylab.xlabel("free energy [kcal/mol]", fontsize=25)
    pylab.ylabel("") #"macro states")
    pylab.savefig(output_file + '.svg')
    pylab.clf()
    pylab.close()
    

class Filter():
    MinHeight = None
    MaxMinima = None
    MaxEnergy = None
    def __init__(self, minh, max, max_e):
        self.MinHeight = minh
        self.MaxMinima = max
        self.MaxEnergy = max_e
        
        
def get_ids_to_color(minima_to_color, structure_list):
    to_color = []
    for m in minima_to_color:
        for s in structure_list:
            if m.Structure == s.Structure:
                to_color.append(s.Index)
                break
    return to_color
 
def read_barriers_structure_map_file(barriers_map_file_path):
    #basin_structure_to_input_line_map = {}
    basin_index_to_input_line_index = {}
    input_line_to_basin_index = {}
    lines = []
    map_regex = "^\s*([\.\(\)]+)\s*(\d+)\s*(\-?\d+\.?\d*)\s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+)\s*(\d+)"
    with open(barriers_map_file_path, 'r') as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        match = re.match(map_regex, line)
        if match:
            basin_representative = match.group(1)
            index_in_energy_sorted_list = match.group(2)
            energy_value_of_input_structure = match.group(3)
            #myms.min, myms.truemin, myms.gradmin, myms.truegradmin
            min_index = match.group(4)
            true_min_index = match.group(5)
            gradient_min_index = match.group(6)
            true_gradient_min_index = match.group(7)
            true_gradient_min_index = int(true_gradient_min_index)
            line_index_in_input_file = match.group(8)
            
            #basin_structure_to_input_line_map[basin_representative] = idx
            if not true_gradient_min_index in basin_index_to_input_line_index:
                basin_index_to_input_line_index[true_gradient_min_index] = set()
            basin_index_to_input_line_index[true_gradient_min_index].add(idx)
            
            if not idx in input_line_to_basin_index:
                input_line_to_basin_index[idx] = set()
            input_line_to_basin_index[idx].add(true_gradient_min_index)
            
    return basin_index_to_input_line_index, input_line_to_basin_index #, basin_structure_to_input_line_map  
    


def barriers_zeugs(saddle_file, filter):
    barriers_minima_representatives = []
    barriers_saddles = []
    barriers_tree_barriers = None
    barriers_saddle_matrix = None
    
    output_file = "bar_tree"
    sequence, barriers_minima_representatives, min_energy, barriers_saddles = read_bar_saddle_file(saddle_file, filter)
    barriers_saddles.sort(key = operator.itemgetter(2))
    print("barriers tree")
    newick_string, barriers_saddles, barriers_tree_barriers, barriers_saddle_matrix = create_newick_tree_string(barriers_saddles, barriers_minima_representatives, filter)
    
    max_saddle_energy = barriers_saddles[-1][2]
    if barriers_saddles[-1][1] == -1: # mfe has not way out to another basin.
        max_saddle_energy = barriers_saddles[-2][2]
    #print(max_saddle_energy, 'max_saddle')
    #print(barriers_saddles)
    ids_to_color = get_ids_to_color(minima_to_color, barriers_minima_representatives)
    #print(ids_to_color)
    #newick_string = "(1:1.0,2:2.0):2,((3:1.0,4:1.0):4);"
    #max_saddle_energy = 6
    #min_energy = 1
    print("max min", max_saddle_energy, min_energy)
    #print_newick_tree(newick_string, output_file, max_saddle_energy, min_energy, ids_to_color)
    
    write_saddle_file(find_path_saddle_list, "barriers_saddles.ssv")
    #for s in barriers_saddles:
    #    print(s)
    return (barriers_minima_representatives, barriers_saddles, barriers_tree_barriers, barriers_saddle_matrix)


def write_matrix_file(numpy_matrix, file_name):
    x_states,y_states = a.shape
    if x_states != y_states:
        print("Error: matrix is not a square!")
    with open(file_name, 'w') as f:
        f.write("# upper triangular matrix")
        for i in range(x_states):
            line = ""
            for j in range(i+1, x_states):
                line += "{:.2f}".format(numpy_matrix[i,j]) + " "
            line += "\n"
            f.write(line)
    return 0

def write_saddle_file(saddles, file_name):
    with open(file_name, 'w') as f:
        for s in saddles:
            index_from = s[0]
            index_to = s[1]
            saddle_energy = s[2]
            line = str(index_from) + " " + str(index_to) + " " + "{:.2f}".format(saddle_energy)
            f.write(line)
    return 0

Max_Threads = 1
       
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read structures and create a barriers tree.')
    parser.add_argument("-f", "--structure_file", type=str, required=False, help="File with RNA secondary structures.")
    parser.add_argument("-s", "--sequence", type=str, required=False, help="Sequence that belongs to the structure file f.")
    parser.add_argument("-b", "--saddle_file", type=str, required=False, help="Barriers output file with saddles.")
    #parser.add_argument("-m", "--minh", type=float, required=False, help="Join minima that are in the range of minimum height.") # not implemented
    parser.add_argument("-x", "--max", type=int, required=False, help="Print the lowest minima.")
    parser.add_argument("-t", "--threads", type=int, required=False, help="Number of threads for parallel find path computations.") 
    parser.add_argument("-e", "--max_energy", type=float, required=False, help="Print everything below this saddle threshold.")
    parser.add_argument("-c", "--color_minima", type=str, required=False, help="File with minima that should be colored (one in each line).")
    parser.add_argument("-z", "--compute_l2_norm_saddle_heights", type=str, required=False, help="Input: barriers file with mapped structures to barriers minima. Needs a structure file for which a barrier tree is constructucted and a barriers file with the reference tree." + \
                        "The basins for the saddles from the artificial tree are mapped into the reference tree then the saddle height in the reference tree is determined."+
                        " The absolute difference of all saddle heights squared is the l2 norm.")
    
    args = parser.parse_args()
    
    sys.setrecursionlimit(int(pow(2,31)-1))
    
    Max_Threads = 1
    if args.threads:
        Max_Threads = args.threads
    
    #filter = Filter(args.minh, args.max, args.max_energy)
    filter = Filter(None, args.max, args.max_energy)
    
    barriers_minima_representatives = []
    barriers_saddles = []
    barriers_tree_barriers = None
    barriers_saddle_matrix = None
    
    gradient_walk_minima_representatives = []
    find_path_saddle_list = []
    barriers_tree_findpath = None
    findpath_saddle_matrix = None
    
    minima_to_color = []
    if args.color_minima:
        minima_to_color = read_pure_structure_list(args.color_minima)
        
    if args.saddle_file:
        barriers_minima_representatives, barriers_saddles, barriers_tree_barriers, barriers_saddle_matrix = barriers_zeugs(args.saddle_file, filter)

    if args.structure_file:
        output_file = "find_path_tree"
        sequence = None
        gradient_walk_minima_representatives = []
        sequence, gradient_walk_minima_representatives, min_energy = read_structure(args.structure_file)

        if filter.MaxMinima != None:
            gradient_walk_minima_representatives.sort(key=lambda x: x.Energy, reverse=False)
            gradient_walk_minima_representatives = gradient_walk_minima_representatives[:filter.MaxMinima]
            
        if sequence == None:
            if args.sequence == None:
                print("Error: file contains no sequence!")
                exit()
         
            else:
                sequence = args.sequence
        
        find_path_saddle_list = connect_structures_find_saddles(sequence, gradient_walk_minima_representatives)
        #sort according to saddle energy
        find_path_saddle_list.sort(key = operator.itemgetter(2))
        
        
        print("findpath tree")
        #create newick string: (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);
        newick_string, find_path_saddle_list, barriers_tree_findpath, findpath_saddle_matrix  = create_newick_tree_string(find_path_saddle_list, gradient_walk_minima_representatives, filter)
        max_saddle_energy = find_path_saddle_list[-1][2]
        if find_path_saddle_list[-1][1] == -1: # mfe has not way out to another basin.
            max_saddle_energy = find_path_saddle_list[-2][2]
        n_minima = len(gradient_walk_minima_representatives)
        ids_to_color = get_ids_to_color(minima_to_color, gradient_walk_minima_representatives)

        #print_newick_tree(newick_string, output_file, max_saddle_energy, min_energy, ids_to_color)
        
        write_saddle_file(find_path_saddle_list, "find_path_saddles.ssv")


    if args.compute_l2_norm_saddle_heights:
        if not args.saddle_file or not args.structure_file:
            print("Error: we need a barriers file for the original tree and a structure file for the find_path tree!")
        barriers_map_file_path = args.compute_l2_norm_saddle_heights
        basin_index_to_input_line_index, input_line_to_basin_index = read_barriers_structure_map_file(barriers_map_file_path)
        
        
        leaf_indices_findpath = barriers_tree_findpath.get_leafs()
        leaf_indices_findpath = list(leaf_indices_findpath)
        
        leaf_indices_barriers = barriers_tree_barriers.get_leafs()
        #for l in leaf_indices_findpath:
        #    input_line_to_basin_index[]
        
        minima_barriers = len(barriers_minima_representatives)
        minima_findpath = len(gradient_walk_minima_representatives)
        minima_barriers_tree = len(leaf_indices_barriers)
        minima_findpath_tree = len(leaf_indices_findpath)
        minima_mapped_from_find_path_to_barriers = len(input_line_to_basin_index.keys())
        minima_mapped_from_barriers_to_findpath = len(basin_index_to_input_line_index.keys())
        minima_mapped_from_find_path_tree_to_barriers_tree = 0
        minima_mapped_from_find_path_tree_to_barriers_tree_unique = set()
        print(leaf_indices_findpath)
        print(input_line_to_basin_index)
        for l in leaf_indices_findpath:
            m = None
            try:
                m = input_line_to_basin_index[l]
            except Exception as e:
                pass # minimum is not in barriers output (maybe maxE filter)
            if m != None:
                for mm in m:
                    if mm in leaf_indices_barriers:
                        minima_mapped_from_find_path_tree_to_barriers_tree += 1
                        minima_mapped_from_find_path_tree_to_barriers_tree_unique.add(mm)
        
        minima_mapped_from_barriers_tree_to_findpath_tree = 0
        for l in leaf_indices_barriers:
            m = None
            try:
                m = basin_index_to_input_line_index[l]
            except Exception as e:
                pass # minimum is not in sampled minima
            if m != None:
                for mm in m:
                    if mm in leaf_indices_findpath:
                        minima_mapped_from_barriers_tree_to_findpath_tree += 1
                  
        all_measures = [minima_barriers, minima_findpath, minima_barriers_tree, minima_findpath_tree, minima_mapped_from_find_path_to_barriers, \
                  minima_mapped_from_barriers_to_findpath, minima_mapped_from_find_path_tree_to_barriers_tree, minima_mapped_from_barriers_tree_to_findpath_tree]
        with open("tree_mapped_indices_statistics.txt", 'w') as f:
            csv_header = "minima_barriers, minima_findpath, minima_barriers_tree, minima_findpath_tree, minima_mapped_from_find_path_to_barriers, " + \
                  "minima_mapped_from_barriers_to_findpath, minima_mapped_from_find_path_tree_to_barriers_tree, minima_mapped_from_barriers_tree_to_findpath_tree"
            f.write(csv_header)
            results_line = ", ".join([ str(x) for x in all_measures])
            f.write(results_line)
            
        with open("percentage_tree_covered.txt", 'w') as f:
            results_line = (len(minima_mapped_from_find_path_tree_to_barriers_tree_unique) / len(leaf_indices_barriers))*100.0
            f.write("{:.2f}".format(results_line))

        #leaf_indices_findpath = barriers_tree_findpath.get_leafs()
        sum_diff = 0
        for i, z_i in enumerate(leaf_indices_findpath):
            set_idx_from = None
            try:
                set_idx_from = input_line_to_basin_index[z_i]
            except Exception as e:
                print("Warning: fp index not in barriers tree", z_i)
                pass

            if set_idx_from != None:
                for z_j in leaf_indices_findpath[i+1:]:
                    barriers_saddle = None
                    findpath_saddle = 0
                    
                    set_idx_to = None
                    try:
                        set_idx_to = input_line_to_basin_index[z_j]
                    except Exception as e:
                        print("Warning: fp index not in barriers tree", z_j)
                        pass
                    
                    if set_idx_to != None:
                        
                        x,y = sorted([z_i,z_j])
                        findpath_saddle = findpath_saddle_matrix[x, y]
                        
                        mapped_saddle_count = 0
                        for m_i in set_idx_from:
                            for m_j in set_idx_to:
                                if m_i != m_j:
                                    a,b = sorted([m_i-1, m_j-1])
                                    if a in leaf_indices_barriers and b in leaf_indices_barriers:
                                        barriers_saddle = barriers_saddle_matrix[a, b]
                                        mapped_saddle_count += 1
                        
                        if mapped_saddle_count > 1:
                            print("Error: fp minima are mapped to several barriers minima! This should not be possible.", z_i, z_i, set_idx_from, set_idx_to)
                        
                        if barriers_saddle != None:
                            saddle_diff = findpath_saddle - barriers_saddle
                            sum_diff += (saddle_diff * saddle_diff)
                    
        l2_norm = sqrt(sum_diff)
        """
        leaf_indices_barriers = barriers_tree_barriers.get_leafs()
        leaf_indices_barriers = list(leaf_indices_barriers)

        #leaf_indices_findpath = barriers_tree_findpath.get_leafs()
        sum_diff = 0
        for i, zb_i in enumerate(leaf_indices_barriers): 
            for zb_j in leaf_indices_barriers[i+1:]:
                barriers_saddle = 0

                x,y = sorted([zb_i,zb_j])
                barriers_saddle= barriers_saddle_matrix[x, y]
                if np.isnan(barriers_saddle):
                    barriers_saddle = 0
                #print("barriers lca:", zb_i, zb_j, barriers_saddle)
                
                set_idx_from = None
                set_idx_to = None
                try:
                    # should be zero based (check the input file!)
                    set_idx_from = basin_index_to_input_line_index[zb_i+1]
                    set_idx_to = basin_index_to_input_line_index[zb_j+1]
                except Exception as e:
                    pass
                print("sets",set_idx_from, set_idx_to, zb_i, zb_j)
                lowest_fp_saddle = sys.maxsize
                if set_idx_from != None and set_idx_to != None:
                    for fp_i in set_idx_from:
                        for fp_j in set_idx_to:
                            if fp_i == fp_j:
                                lowest_fp_saddle = 0
                                break
                            else:
                                #fp_saddle = barriers_tree_findpath.find_lca_saddle_energy(fp_i, fp_j)
                                fp_saddle = 0
                                try:
                                    fp_saddle = findpath_saddle_matrix[fp_i, fp_j]
                                except Exception as e:
                                    pass
                                print("fp lca:", fp_i, fp_j, fp_saddle)
                                if np.isnan(fp_saddle):
                                    fp_saddle = 0
                                if fp_saddle < lowest_fp_saddle:
                                    lowest_fp_saddle = fp_saddle
                else:
                    lowest_fp_saddle = 0 # TODO: handle the case of missing mapped values.
                
                if lowest_fp_saddle == sys.maxsize:
                    print("Error: saddle not found!", b_i, b_j)
                
                # TODO: find proper value if barriers saddle is not in the tree with the mfe or sampled minimum is not in the tree with the mfe.
                if barriers_saddle == None:
                    barriers_saddle = 0
                if lowest_fp_saddle == None:
                    lowest_fp_saddle = 0
                    
                saddle_diff = barriers_saddle - lowest_fp_saddle
                sum_diff += (saddle_diff * saddle_diff)
                
        l2_norm = sqrt(sum_diff)
        """
        
        with open("l2_norm_direct_neighbor_basins.txt", 'w') as f:
            f.write(str(l2_norm))
            
            
            
            
            
            
            
            
            
            
        
        
   

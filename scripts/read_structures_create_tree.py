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
                print(index, structure, energy, structure_basin_index,saddle_energy_distance, saddle_structure)
                if energy < min_energy:
                    min_energy = energy
                try:
                    index = int(index)
                except Exception as e:
                    new_index += 1
                    index = new_index
                structure_object = StructureEnergy(index, structure, energy)
                structure_list.append(structure_object)
                print(index, structure, energy, saddle_energy_kcal)
                saddle_list.append((index-1, structure_basin_index-1, saddle_energy_kcal))
            else:
                match = re.match("\s*([ACGUTNIP]+)", line)
                if match:
                    sequence = match.group(1)
    #print(saddle_list)
    return sequence, structure_list, min_energy, saddle_list


def connect_structures_find_saddles(sequence, structure_list):
    pairs = {}
    fc = RNA.fold_compound(sequence)
    for i, se_1 in enumerate(structure_list):
        for j in range(i+1, len(structure_list)):
            se_2 = structure_list[j]
            saddle_energy_dcal = fc.path_findpath_saddle(se_1.Structure, se_2.Structure)
            #saddle_energy_dcal_2 = fc.path_findpath_saddle(se_2.Structure, se_1.Structure)
            #saddle_energy_dcal = min(saddle_energy_dcal, saddle_energy_dcal_2)
            saddle_energy_kcal = saddle_energy_dcal / 100.0
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
    Parent = None
    Branch_Length = 0.0
    Saddle_Energy = None
    Energy = None
    Index = None
    def __init__(self, index, energy, saddle_energy, parent = None, children = None):
        self.Index = index
        self.Energy = energy
        self.Saddle_Energy = saddle_energy
        if saddle_energy != None and energy != None:
            self.Branch_Length = self.Saddle_Energy - self.Energy
        self.Parent = parent
        if type(children) == type([]):
            self.Children = children
    
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
            

def create_barrier_tree(minimal_saddle_list, structure_list, filter, do_merge_ignore_connectivity = False):
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
        print(id_from, id_to)
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
                    print("update lengths", id_from, id_to)
                    # update lengths only if they are not in the same cluster
                    clusters[c_id_from].extend_branch_length(saddle_difference)
                    clusters[c_id_to].extend_branch_length(saddle_difference)
                    # and merge them!
                    
                    parent = Node(None, None, saddle_energy, None, children = [clusters[c_id_from], clusters[c_id_to]])
                    clusters.pop(max(c_id_from, c_id_to))
                    clusters.pop(min(c_id_from, c_id_to))
                    clusters.append(parent)
                else:
                    # update lengths in same cluster
                    clusters[c_id_from].extend_branch_length(saddle_difference)
                    
                has_updated_length_for_ids.append(id_pair)
            continue
        
        if not id_from in clustered_ids and not id_to in clustered_ids:
            print("new cluster", id_from, id_to)
            # create new cluster with both ids
            energy_from = structure_list[id_from].Energy
            energy_to = structure_list[id_to].Energy
            node_from = Node(id_from, energy_from, saddle_energy, parent = None, children = None)
            node_to = Node(id_to, energy_to, saddle_energy, parent = None, children = None)
            
            index = None
            energy = None
            #saddle_energy = None
            parent = Node(index, energy, saddle_energy, parent = None, children = [node_from, node_to])
            node_from.Parent = parent
            node_to.Parent = parent
            
            if do_merge_ignore_connectivity:
                # merge previous clusters and add new cluster
                if len(clusters) > 1:
                    index_merge = None
                    energy_merge = None
                    #saddle_energy_merge = None
                    copied_children = [ x for x in clusters ]
                    del clusters[:]
                    merged_clusters = Node(index_merge, energy_merge, saddle_energy, parent = None, children = copied_children)
                    merged_clusters.Branch_Length = saddle_difference
                    clusters.append(merged_clusters)
            
            for c in clusters:
                #c.extend_branch_length(saddle_difference)
                diff = saddle_energy - c.Saddle_Energy
                c.Branch_Length += diff
                c.Saddle_Energy = saddle_energy
            
            clusters.append(parent)
                
            clustered_ids.add(id_from)
            clustered_ids.add(id_to)
            continue
        
        if id_from in clustered_ids and not id_to in clustered_ids:
            print("add cluster", id_to)
            
            node_from = None
            if not do_merge_ignore_connectivity:
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
            else:
                # create node that connects all clusters (does not matter if they contain id_from)
                index_merge = None
                energy_merge = None
                #saddle_energy_merge = None
                copied_children = [ x for x in clusters ]
                del clusters[:]
                node_from = Node(index_merge, energy_merge, saddle_energy, parent = None, children = copied_children)
                node_from.Branch_Length = saddle_difference
            
            energy_to = structure_list[id_to].Energy
            node_to = Node(id_to, energy_to, saddle_energy, parent = None, children = None)
            
            index = None
            energy = None
            #saddle_energy = None
            parent = Node(index, energy, saddle_energy, parent = None, children = [node_from, node_to])
            
            node_from.Parent = parent
            node_to.Parent = parent
            
            for c in clusters:
                #c.extend_branch_length(saddle_difference)
                diff = saddle_energy - c.Saddle_Energy
                c.Branch_Length = +diff
                c.Saddle_Energy = saddle_energy
                
            clusters.append(parent)
            
            clustered_ids.add(id_to)
            continue
            
        if not id_from in clustered_ids and id_to in clustered_ids:
            print("add cluster", id_from)
            
            node_to = None
            if not do_merge_ignore_connectivity:
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
            else:
                # create node that connects all clusters (does not matter if they contain id_to)
                index_merge = None
                energy_merge = None
                #saddle_energy_merge = None
                copied_children = [ x for x in clusters ]
                del clusters[:]
                node_to = Node(index_merge, energy_merge, saddle_energy, parent = None, children = copied_children)
                node_to.Branch_Length = saddle_difference

            energy_from = structure_list[id_from].Energy
            node_from = Node(id_from, energy_from, saddle_energy, parent = None, children = None)
            
            index = None
            energy = None
            #saddle_energy = None
            parent = Node(index, energy, saddle_energy, parent = None, children = [node_from, node_to])
            
            node_from.Parent = parent
            node_to.Parent = parent
            
            for c in clusters:
                #c.extend_branch_length(saddle_difference)
                diff = saddle_energy - c.Saddle_Energy
                c.Branch_Length += diff
                c.Saddle_Energy = saddle_energy
            
            clusters.append(parent)
            
            clustered_ids.add(id_from)
            continue
    
    return clusters


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
    
  
def create_newick_tree_string(minimal_saddle_list, structure_list, filter, do_merge_ignore_connectivity=False):
    # create barriers tree (could be one or many if it is unconnected)
    clusters = create_barrier_tree(minimal_saddle_list, structure_list, filter, do_merge_ignore_connectivity)
    
    if len(clusters) <= 0:
        print("Error: no tree was created!")
    if len(clusters) > 1:
        print("Error: the landscape is not connected! We have a forest!")
    print(len(clusters),"cl")
    barriers_tree = None
    for i in range(len(clusters)):
        if clusters[i].contains_node_index(0):
            # if it contains the mfe basin (with index 0)
            barriers_tree = clusters.pop(i)
            break

    print(len(clusters),"cl")
    newick_string = ""
    if barriers_tree != None:
        newick_string = newick_string_builder(barriers_tree) + ";"
        print("tree to plot:")
        print(newick_string)
    else:
        print("Error: the mfe structure is not connected to the tree!")
    
    if len(clusters) > 0:
        print("other trees")
        for c in clusters:
            other_tree = newick_string_builder(c) + ";"
            print(other_tree)
            if newick_string == "":
                newick_string = other_tree
    print(len(clusters),"cl")

    return newick_string, minimal_saddle_list


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
                
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read structures and create a barriers tree.')
    parser.add_argument("-f", "--structure_file", type=str, required=False, help="File with RNA secondary structures.")
    parser.add_argument("-b", "--saddle_file", type=str, required=False, help="Barriers output file with saddles.")
    parser.add_argument("-m", "--minh", type=float, required=False, help="Join minima that are in the range of minimum height.")
    parser.add_argument("-x", "--max", type=int, required=False, help="Print the lowest minima.") 
    parser.add_argument("-e", "--max_energy", type=float, required=False, help="Print everything below this saddle threshold.")
    parser.add_argument("-c", "--color_minima", type=str, required=False, help="File with minima that should be colored (one in each line).")
    args = parser.parse_args()

    filter = Filter(args.minh, args.max, args.max_energy)
    
    minima_to_color = []
    if args.color_minima:
        sequence, minima_to_color, min_energy = read_structure(args.color_minima)
    if args.saddle_file:
        output_file = "bar_tree"
        sequence, structure_list, min_energy, saddle_list = read_bar_saddle_file(args.saddle_file, filter)
        saddle_list.sort(key = operator.itemgetter(2))
        newick_string, saddle_list = create_newick_tree_string(saddle_list, structure_list, filter)
        
        max_saddle_energy = saddle_list[-1][2]
        if saddle_list[-1][1] == -1: # mfe has not way out to another basin.
            max_saddle_energy = saddle_list[-2][2]
        print(max_saddle_energy, 'max_saddle')
        print(saddle_list)
        ids_to_color = get_ids_to_color(minima_to_color, structure_list)
        print(ids_to_color)
        #newick_string = "(1:1.0,2:2.0):2,((3:1.0,4:1.0):4);"
        #max_saddle_energy = 6
        #min_energy = 1
        print_newick_tree(newick_string, output_file, max_saddle_energy, min_energy, ids_to_color)
        for s in saddle_list:
            print(s)

    if args.structure_file:
        output_file = "find_path_tree"
        sequence = None
        structure_list = []
        sequence, structure_list, min_energy = read_structure(args.structure_file)
        
        if filter.MaxMinima != None:
            structure_list.sort(key=lambda x: x.Energy, reverse=False)
            structure_list = structure_list[:filter.MaxMinima]
            
        if sequence == None:
            print("Error: file contains no sequence!")
            exit()
        
        minimal_saddle_list = connect_structures_find_saddles(sequence, structure_list)
        #sort according to saddle energy
        minimal_saddle_list.sort(key = operator.itemgetter(2))

        #create newick string: (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);
        newick_string, minimal_saddle_list  = create_newick_tree_string(minimal_saddle_list, structure_list, filter, do_merge_ignore_connectivity=False)
        max_saddle_energy = minimal_saddle_list[-1][2]
        if minimal_saddle_list[-1][1] == -1: # mfe has not way out to another basin.
            max_saddle_energy = minimal_saddle_list[-2][2]
        n_minima = len(structure_list)
        ids_to_color = get_ids_to_color(minima_to_color, structure_list)

        print_newick_tree(newick_string, output_file, max_saddle_energy, min_energy, ids_to_color)

        
        
        
        
        
        
        
        
   

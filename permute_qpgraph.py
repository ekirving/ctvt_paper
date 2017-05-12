#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Branch and bound

# from Bio import Phylo
import xml.etree.ElementTree as ET
import copy, hashlib, re

# import my custom modules
from pipeline_utils import *

MAX_OUTLIER_THRESHOLD = 0

parfile = 'permute/permute.par'

root = 'R'
out = 'Out'

# nodes = ['A', 'B']  # 2
# nodes = ['A', 'B', 'C']  # 3
nodes = ['A', 'B', 'C', 'X']  # 4
# nodes = ['A', 'B', 'C', 'X', 'Q']  # 5
# nodes = ['A', 'B', 'C', 'X', 'Q', 'D']  # 6
# nodes = ['A', 'B', 'C', 'X', 'Q', 'D', 'S']  # 7


def add_node(root_tree, new_node):

    # get all the leaf nodes
    nodes = root_tree.findall('.//*')

    for node in nodes:
        # skip the outgroup and non-leaf nodes
        if node.tag == out or len(node) > 0:
            continue

        # clone the current tree
        new_tree = copy.deepcopy(root_tree)

        # get the target node in the new tree
        target_node = new_tree.find('.//' + node.tag)

        # get the parent of the target
        parent_node = new_tree.find('.//' + node.tag + '/..')

        # does the target node have a sibling
        if len(parent_node) > 1:
            parent_node.remove(target_node)

            # add an intermediate node, to act as the parent for the new node
            parent_node = ET.SubElement(parent_node, new_label())

            # move the target node
            ET.SubElement(parent_node, target_node.tag)

        # add the new node as a sibling to the target
        ET.SubElement(parent_node, new_node)

        yield new_tree


def build_tree(root_tree, unplaced):

    for new_node in unplaced:

        # adding each new node generates (n-1)*2 new trees, where n is the number of nodes already in the tree
        new_trees = add_node(root_tree, new_node)

        # get the remaining nodes
        remaining = list(unplaced)
        remaining.remove(new_node)

        for new_tree in new_trees:

            # convert the tree to qpGraph format
            graph = convert_tree(new_tree)

            # get a unique names for the output files
            graph_name = hash_text(graph)
            grp_file = 'permute/graphs/{name}.graph'.format(name=graph_name)
            dot_file = 'permute/graphs/{name}.dot'.format(name=graph_name)
            log_file = 'permute/graphs/{name}.log'.format(name=graph_name)

            # save the graph file
            with open(grp_file, 'w') as fout:
                fout.write(graph)

            # run qpGraph
            log = run_cmd(["qpGraph",
                           "-p", parfile,
                           "-g", grp_file,
                           "-d", dot_file])

            # save the log file
            with open(log_file, 'w') as fout:
                fout.write(log)

            # parse the log and extract the outliers
            outliers, worst_fstat = extract_outliers(log)

            print "\t".join(worst_fstat)

            # for each new tree that passes threshold, lets add the remaining nodes
            if len(outliers) <= MAX_OUTLIER_THRESHOLD:
                build_tree(new_tree, remaining)


def extract_outliers(log):

    outliers = []
    read_log = False
    worst_fstat = ['No outliers!']

    for line in log:
        if 'outliers' in line:
            read_log = True
            continue
        elif 'worst f-stat' in line:
            worst_fstat = line.split()
            read_log = False
            continue

        if read_log and len(line.strip()) > 0:
            # save all the outliers
            outliers.append(line.split())

    return outliers, worst_fstat


def convert_tree(root_tree):

    graph = "root\t{root}\n".format(root=root)

    # get all the leaf nodes
    nodes = root_tree.findall('.//*')

    for node in nodes:
        # skip non-leaf nodes
        if len(node) == 0:
            graph += "label\t{node}\t{node}\n".format(node=node.tag)

    graph += convert_node(root_tree.getroot())

    return graph


def convert_node(parent_node):
    graph = ""

    for child_node in parent_node:

        child_code = hash_text(child_node.tag)

        graph += "edge\t{code}\t{parent}\t{inter}\n".format(code=child_code, parent=parent_node.tag, inter=child_node.tag)

        # leaf nodes
        if len(child_node) > 0:
            # now convert the children
            graph += convert_node(child_node)

    return graph


def hash_text(text, len=7):
    return hashlib.sha1(text).hexdigest()[0:len]


# TODO tidy this up
def new_label():
    """
    Return a new label for a node
    """
    new_label.n += 1
    return 'n{}'.format(new_label.n)
new_label.n = 0

unplaced = list(nodes)

# lets get started
for node in nodes:

    # setup the simplest 2-node tree
    root_node = ET.Element(root)
    root_tree = ET.ElementTree(root_node)

    ET.SubElement(root_node, out)
    ET.SubElement(root_node, node)

    # get the unplaced nodes
    unplaced.remove(node)

    build_tree(root_tree, unplaced)
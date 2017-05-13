#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Branch and Bound Algorithm /
# Randomised Stepwise Addition Order Algorithm

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


def permute_tree(root_tree, new_node):
    """
    Permute all possible new trees, by adding the new node to all external branches and all possible pairs of external
    branches.
     
    :param root_tree: The tree to permute 
    :param new_node: The node to add
    :return: One of the permuted trees
    """

    # get all the leaf nodes
    nodes = root_tree.findall('.//*')

    # add the new node to every leaf node in the tree
    for target_node in list(nodes):

        # skip the outgroup and non-leaf nodes
        if target_node.tag == out or len(target_node) > 0:
            nodes.remove(target_node)
            continue

        # clone the current tree and add the node
        new_tree = copy.deepcopy(root_tree)
        insert_node(new_tree, target_node, new_node)

        yield new_tree

    if len(nodes) > 1:
        # now permute all the two parent admixture possibilities
        pairs = list(itertools.combinations(nodes, 2))

        for target1, target2 in pairs:

            # clone the current tree
            new_tree = copy.deepcopy(root_tree)

            # add the new node as the child of both targets
            insert_node(new_tree, target1, new_node)
            insert_node(new_tree, target2, new_node)

            yield new_tree


def insert_node(new_tree, target_node, new_node):
    """
    Helper function to add a new node to the branch leading to the target node.
    
    :param new_tree: The tree to add the node to 
    :param target_node: The nodes attached to the target branch
    :param new_node: The new node to add
    :return: The new tree
    """
    # get the target node in the new tree
    target_node = new_tree.find('.//' + target_node.tag)

    # get the parent of the target
    parent_node = new_tree.find('.//' + target_node.tag + '/..')

    # does the target node have a sibling
    if len(parent_node) > 1:
        parent_node.remove(target_node)

        # add an intermediate node, to act as the parent for the new node
        parent_node = ET.SubElement(parent_node, new_label())

        # re add the target node
        ET.SubElement(parent_node, target_node.tag)

    # add the new node as a sibling to the target
    ET.SubElement(parent_node, new_node)


def recurse_tree(root_tree, unplaced):

    for new_node in unplaced:

        # permute all possible new trees and graphs
        new_trees = permute_tree(root_tree, new_node)

        # get the remaining nodes
        remaining = list(unplaced)
        remaining.remove(new_node)

        for new_tree in new_trees:

            # run qpGraph to test the model
            outliers = run_qpgraph(new_tree)

            # for each new tree that passes threshold, let's add the remaining nodes
            if len(outliers) <= MAX_OUTLIER_THRESHOLD:
                if len(remaining) == 0:
                    print "SUCCESS: Placed all nodes on a graph without outliers!"
                else:
                    recurse_tree(new_tree, remaining)


def run_qpgraph(new_tree):

    # export newick tree
    newick = export_newick_tree(new_tree.getroot()).ljust(len(nodes) * 7)

    # convert the tree to qpGraph format
    graph = export_qpgraph(new_tree)

    # get a unique names for the output files
    graph_name = hash_text(graph)
    grp_file = 'permute/graphs/{name}.graph'.format(name=graph_name)
    dot_file = 'permute/graphs/{name}.dot'.format(name=graph_name)
    log_file = 'permute/graphs/{name}.log'.format(name=graph_name)
    pdf_file = 'permute/graphs/{name}.pdf'.format(name=graph_name)

    # save the graph file
    with open(grp_file, 'w') as fout:
        fout.write(graph)

    # run qpGraph
    log = run_cmd(["qpGraph", "-p", parfile, "-g", grp_file, "-d", dot_file], verbose=False)

    # make the PDF
    pdf = run_cmd(["dot", "-Tpdf", dot_file], verbose=False)

    # save the ps file
    with open(pdf_file, 'w') as fout:
        fout.write(pdf)

    # save the log file
    with open(log_file, 'w') as fout:
        fout.write(log)

    # parse the log and extract the outliers
    outliers, worst_fstat = extract_outliers(log.splitlines())

    # print some summary stats
    print "{name}\ttree={tree}\toutliers={out}\tworst={worst}".format(name=graph_name, tree=newick,
                                                                      out=len(outliers), worst=worst_fstat[-1])

    return outliers


def extract_outliers(log):

    outliers = []
    read_log = False
    worst_fstat = []

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


def export_qpgraph(root_tree):

    graph = "root\t{root}\n".format(root=root)

    # get all the nodes
    nodes = root_tree.findall('.//*')

    for node in nodes:
        # skip non-leaf nodes
        if len(node) == 0:
            graph += "label\t{node}\t{node}\n".format(node=node.tag)

    graph += export_qpgraph_node(root_tree)

    return graph


def export_qpgraph_node(root_tree, parent_node=None):
    graph = ""

    if parent_node is None:
        parent_node = root_tree.getroot()

    for child_node in parent_node:

        # is this an admixture node or a normal node
        matches = root_tree.findall('.//' + child_node.tag + '/..')

        if len(matches) > 1:
            # admixture branch
            new_node = new_label()
            parent1, parent2 = matches
            graph += "admix\t{child}\t{parent1}\t{parent2}\t50\t50\n".format(parent1=parent1.tag,
                                                                             parent2=parent2.tag,
                                                                             child=new_node)

            # remove both nodes so we don't export them twice
            parent1.remove(parent1.find(child_node.tag))
            parent2.remove(parent2.find(child_node.tag))

            parent_node.tag = new_node

        # regular branch
        code = hash_text(child_node.tag)
        graph += "edge\t{code}\t{parent}\t{child}\n".format(code=code, parent=parent_node.tag, child=child_node.tag)

        # leaf nodes
        if len(child_node) > 0:
            # now convert the children
            graph += export_qpgraph_node(root_tree, child_node)

    return graph


def hash_text(text, len=7):
    return hashlib.sha1(text).hexdigest()[0:len]


def new_label():
    """
    Return a new label for a node
    """
    new_label.n += 1
    return 'n{}'.format(new_label.n)
new_label.n = 0


def export_newick_tree(parent_node):
    """
    Convert an XML tree into Newick format
    """
    if len(parent_node) == 0:
        return parent_node.tag
    else:
        children = [export_newick_tree(child_node) for child_node in parent_node]
        return '(' + ','.join(children) + ')%s' % parent_node.tag


unplaced = list(nodes)

# lets get started
for node in nodes:

    # initialise all of the simplest 2-node trees
    root_node = ET.Element(root)
    root_tree = ET.ElementTree(root_node)

    ET.SubElement(root_node, out)
    ET.SubElement(root_node, node)

    # get the unplaced nodes
    unplaced.remove(node)

    recurse_tree(root_tree, unplaced)
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Build all possible trees and graphs using a Randomised Stepwise Addition Order Algorithm

import xml.etree.ElementTree as ElemTree
import re

import itertools
from multiprocessing import Pool

# import the custom modules
from pipeline_utils import *

MULTITHREAD_SEARCH = True

MAX_OUTLIER_THRESHOLD = 0

OUTPUT_FOLDER = 'permute/graphs/'

# parfile = 'permute/permute.par'
# ROOT = 'R'
# OUT = 'Out'
# NODES = ['A', 'B', 'C', 'X']

parfile = 'permute/merged_v2_hq2_nomex_ctvt.par'
root = 'R'
OUT = 'WAM'
NODES = ['DEU', 'DCH', 'DPC', 'CTVT', 'DHU', 'DGL', 'DMA']


def permute_tree(root_tree, new_node):
    """
    Permute all possible new trees, by adding the new node to all external branches and all possible pairs of external
    branches.
    """
    # get all the leaf nodes
    target_nodes = list(root_tree.findall('.//*'))

    # add the new node to every leaf node in the tree
    for target_node in list(target_nodes):

        # skip the outgroup and non-leaf nodes
        if target_node.tag == OUT or len(target_node) > 0:
            target_nodes.remove(target_node)
            continue

        # clone the current tree and add the node
        new_tree = copy.deepcopy(root_tree)
        insert_node(new_tree, target_node, new_node)

        yield new_tree

    # TODO change code to only test admix brances if we can't place a node properly
    if len(target_nodes) > 1:
        # now permute all the two parent admixture possibilities
        pairs = list(itertools.combinations(target_nodes, 2))

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
    """
    # get the target node in the new tree
    target_node = new_tree.find('.//' + target_node.tag)

    # get the parent of the target
    parent_node = new_tree.find('.//' + target_node.tag + '/..')

    # does the target node have a sibling
    if len(parent_node) > 1:
        parent_node.remove(target_node)

        # add an intermediate node, to act as the parent for the new node
        parent_node = ElemTree.SubElement(parent_node, new_label())

        # re add the target node
        ElemTree.SubElement(parent_node, target_node.tag)

    # add the new node as a sibling to the target
    ElemTree.SubElement(parent_node, new_node)


def recurse_tree(root_tree, unplaced):
    """
    Recursively add all the unplaced nodes to the given tree 
    """
    for new_node in unplaced:

        # permute all possible new trees and graphs
        new_trees = permute_tree(root_tree, new_node)

        # get the remaining nodes
        remaining = list(unplaced)
        remaining.remove(new_node)

        for new_tree in new_trees:

            # run qpGraph to test the model
            outliers, num_nodes = run_qpgraph(new_tree)

            # for each new tree that passes threshold, let's add the remaining nodes
            if len(outliers) <= MAX_OUTLIER_THRESHOLD:
                if num_nodes == len(NODES)+1:
                    print "\tSUCCESS: Placed all nodes on a graph without outliers!"
                else:
                    recurse_tree(new_tree, remaining)


def run_qpgraph(new_tree):
    """
    Run qpGraph on the given tree 
    """
    # convert the tree to newick format
    newick = export_newick_tree(new_tree.getroot()).ljust(90)

    # get unique names for the output files
    graph_name = hash_text(newick)
    grp_file = OUTPUT_FOLDER + '{name}.graph'.format(name=graph_name)
    dot_file = OUTPUT_FOLDER + '{name}.dot'.format(name=graph_name)
    log_file = OUTPUT_FOLDER + '{name}.log'.format(name=graph_name)

    cached = False

    try:
        # if the log file exists then we've run the analysis already
        with open(log_file, 'r') as fin:
            log = fin.read()

        cached = True

    except IOError:
        # convert the tree to qpGraph format
        graph = export_qpgraph(new_tree)

        # save the graph file
        with open(grp_file, 'w') as fout:
            fout.write(graph)

        # run qpGraph
        log = run_cmd(["qpGraph", "-p", parfile, "-g", grp_file, "-d", dot_file], verbose=False)

        # save the log file
        with open(log_file, 'w') as fout:
            fout.write(log)

    # parse the log and extract the outliers
    outliers, worst_fstat = extract_outliers(log.splitlines())

    # count the leaf nodes
    leaf_nodes = [node.tag for node in new_tree.findall('.//*') if len(node) == 0]
    num_nodes = len(set(leaf_nodes))
    num_admix = (len(leaf_nodes)-num_nodes)
    num_outliers = len(outliers)

    # embed some useful info in the PDF name
    pdf_file = OUTPUT_FOLDER + '{name}-n{nodes}-o{out}-a{admix}.pdf'.format(name=graph_name,
                                                                           nodes=num_nodes,
                                                                           out=num_outliers,
                                                                           admix=num_admix)

    # TODO and not cached
    if num_outliers <= MAX_OUTLIER_THRESHOLD:
        # generate the PDF
        pdf = run_cmd(["dot", "-Tpdf", dot_file], verbose=False)

        # save the pdf file
        with open(pdf_file, 'w') as fout:
            fout.write(pdf)

    if not cached:
        # print some summary stats
        print "{name}\t{tree}\tnodes={nodes}\tadmix={admix}\toutliers={out}\tworst={worst}".format(
            name=graph_name, tree=newick, nodes=num_nodes, admix=num_admix, out=len(outliers), worst=worst_fstat[-1])

    return outliers, num_nodes


def extract_outliers(log):
    """
    Parse the log file and extract the outliers
    """
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
    """
    Convert the ElementTree into qpGraph format
    """
    # clone the tree because this process is destructive
    local_tree = copy.deepcopy(root_tree)

    graph = "root\t{root}\n".format(root=ROOT)

    # get all the nodes
    nodes = local_tree.findall('.//*')

    for node in nodes:
        # skip non-leaf nodes
        if len(node) == 0:
            graph += "label\t{node}\t{node}\n".format(node=node.tag)

    # build the list of edges
    graph += export_qpgraph_node(local_tree)

    return graph


def export_qpgraph_node(root_tree, parent_node=None):
    """
    Recursively export all the edges in the graph
    """
    graph = ""

    if parent_node is None:
        parent_node = root_tree.getroot()

    for child_node in list(parent_node):

        if child_node.get('printed') == '1':
            continue

        # is this an admixture node or a normal node
        matches = root_tree.findall('.//' + child_node.tag + '/..')

        parent_name = parent_node.tag

        if len(matches) > 1:
            # admixture branch
            new_node = new_label()
            parent1, parent2 = matches
            graph += "admix\t{child}\t{parent1}\t{parent2}\t50\t50\n".format(parent1=parent1.tag,
                                                                             parent2=parent2.tag,
                                                                             child=new_node)

            # flag both nodes so we don't export them twice
            parent1.find(child_node.tag).set('printed', '1')
            parent2.find(child_node.tag).set('printed', '1')

            parent_name = new_node

        # regular branch
        code = hash_text(child_node.tag)
        graph += "edge\t{code}\t{parent}\t{child}\n".format(code=code, parent=parent_name, child=child_node.tag)

        # leaf nodes
        if len(child_node) > 0:
            # now convert the children
            graph += export_qpgraph_node(root_tree, child_node)

    return graph


def hash_text(text, length=7):
    """
    Generate a unique key by hashing a string
    """
    return hashlib.sha1(text).hexdigest()[0:length]


def new_label():
    """
    Return a new label for a node
    """
    new_label.n += 1
    return 'n{}'.format(new_label.n)
new_label.n = 0


def export_newick_tree(parent_node):
    """
    Convert an ElementTree tree into Newick format
    """
    if len(parent_node) == 0:
        return parent_node.tag
    else:
        children = [(child_node.tag ,export_newick_tree(child_node)) for child_node in parent_node]
        children.sort()
        tag_name = '' if re.match('n[0-9]+|R', parent_node.tag) else parent_node.tag
        return '(' + ','.join(node for tag, node in children) + ')%s' % tag_name


def run_analysis(all_nodes):
    """
    Build and test all possible trees and graphs
    """
    data = []
    unplaced = list(all_nodes)

    # make a nested list of all the nodes to place
    for node in all_nodes:
        # e.g. (A,[B,C,D]), (B,[C,D]), (C,[D])
        unplaced.remove(node)
        data.append(list(unplaced))

    if MULTITHREAD_SEARCH:
        with Pool(MAX_CPU_CORES) as pool:
            # initialise all of the simplest 2-node trees
            pool.map(initialise_tree, itertools.izip(all_nodes, data))
    else:
        # initialise trees without multi-threading
        for args in itertools.izip(all_nodes, data):
            initialise_tree(args)


def initialise_tree(args):
    """
    Setup a simple 2-node tree, then recursively add all the unplaced nodes.
    """
    node, unplaced = args

    root_node = ElemTree.Element(ROOT)
    root_tree = ElemTree.ElementTree(root_node)

    ElemTree.SubElement(root_node, OUT)
    ElemTree.SubElement(root_node, node)

    recurse_tree(root_tree, unplaced)

# perform the analysis
run_analysis(NODES)

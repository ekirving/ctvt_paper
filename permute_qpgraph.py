#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Build all possible trees and graphs using a Randomised Stepwise Addition Order Algorithm

import xml.etree.ElementTree as ElemTree
import re

from multiprocessing import Pool

# import the custom modules
from pipeline_utils import *

# MULTITHREAD_SEARCH = True
MULTITHREAD_SEARCH = False

MAX_OUTLIER_THRESHOLD = 0
PROBLEM_NODES = []

PAR_FILE = 'permute/simulated.par'
OUTPUT_FOLDER = 'permute/simulated/'
ROOT = 'R'
OUT = 'Out'
NODES = ['A', 'B', 'X', 'C']


# PAR_FILE = 'permute/merged_v2_hq2_nomex_ctvt.par'
# OUTPUT_FOLDER = 'permute/graphs/'
# ROOT = 'R'
# OUT = 'WAM'
# NODES = ['DEU', 'DCH', 'DPC', 'CTVT', 'DHU', 'DGL', 'DMA']


def recurse_tree(root_tree, new_tag, remaining, depth=0):
    """
    Permute all possible new trees, by adding the new node to all branches.

    If no resulting tree passes the outlier threshold then try adding the node to all possible pairs of branches.
    """
    new_trees = []

    # get all the nodes in the tree (skip the outgroup)
    target_nodes = [node for node in root_tree.findall('.//*') if node.tag != OUT]

    # add the new node to every branch in the tree
    for target_node in target_nodes:

        # clone the current tree and add the new node
        new_tree = copy.deepcopy(root_tree)
        insert_node(new_tree, target_node, new_tag)
        new_trees.append(new_tree)

    # test all the trees
    results = test_trees(new_trees, depth)

    # process the results
    node_placed = check_results(results, remaining, depth)

    # test all the admixture possibilities
    if not node_placed:

        admix_trees = []

        # permute all the two parent admixture possibilities
        pairs = list(itertools.combinations(target_nodes, 2))

        print "\nPairs to check admixture..."

        for target1, target2 in pairs:
            # TODO fixme
            print target1.tag + " " + target2.tag

            # TODO handle targets with matching tag names
            # replace the two targets with an intermediary node
            # add one target as a child to the intermediate node
            # add the new node as a sibling to the target node

            # clone the current tree
            new_tree = copy.deepcopy(root_tree)

            # add the new node as the child of both targets
            insert_node(new_tree, target1, new_tag)
            insert_node(new_tree, target2, new_tag)

            admix_trees.append(new_tree)

        print "\n"

        # test all the admixture trees
        results = test_trees(admix_trees, depth)

        # process the results
        node_placed = check_results(results, remaining, depth)

    # TODO if still not placed then test inserting into a branch
    # TODO should not admix from a branch

    if not node_placed:

        # we could not place the node via either method :(
        if new_tag not in PROBLEM_NODES and remaining:
            print "\tWARNING: Unable to place node '%s' at this time." % new_tag

            PROBLEM_NODES.append(new_tag)

            # add the problem node to end of the list, as we may be able to add it later on
            remaining.append(new_tag)

            # try and add the other nodes
            recurse_tree(root_tree, remaining[0], remaining[1:], depth)

        else:
            print "\tERROR: Cannot place node '%s' in the graph." % new_tag


def test_trees(new_trees, depth):
    """
    Run qpGraph on a list of trees
    """
    if MULTITHREAD_SEARCH:
        # we need to buffer the results to use multi-threading
        pool = Pool(MAX_CPU_CORES)
        results = pool.map(run_qpgraph, itertools.izip(new_trees, itertools.repeat(depth)))
    else:
        # test the trees without multi-threading
        results = []
        for new_tree in new_trees:
            result = run_qpgraph((new_tree, depth))
            results.append(result)

    return results


def check_results(results, remaining, depth):
    """
    Check the results from qpGraph
    """
    # were we able to place the new node
    placed_node = False

    for new_tree, outliers in results:

        # did our new trees pass the threshold
        if len(outliers) <= MAX_OUTLIER_THRESHOLD:

            # recursively add any remaining nodes
            if remaining:
                recurse_tree(new_tree, remaining[0], remaining[1:], depth + 1)
            else:
                print "\tSUCCESS: Placed all nodes on a graph without outliers!"

            # we successfully placed the new node!
            placed_node = True

    return placed_node


def insert_node(new_tree, target_node, new_tag):
    """
    Helper function to add a new node on the branch leading to the target node.
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
        # ElemTree.SubElement(parent_node, target_node.tag)
        parent_node.append(target_node)

    # add the new node as a sibling to the target
    ElemTree.SubElement(parent_node, new_tag)


def run_qpgraph(args):
    """
    Run qpGraph on the given tree 
    """
    # extract the tuple of arguments
    new_tree, depth = args

    # convert the tree to newick format
    newick = export_newick_tree(new_tree.getroot())

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
        log = run_cmd(["qpGraph", "-p", PAR_FILE, "-g", grp_file, "-d", dot_file], verbose=False)

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
    pdf_file = OUTPUT_FOLDER + 'qpgraph-n{nodes}-o{out}-a{admix}-{name}.pdf'.format(nodes=num_nodes,
                                                                                    out=num_outliers,
                                                                                    admix=num_admix,
                                                                                    name=graph_name)

    # TODO fixme
    # if num_outliers <= MAX_OUTLIER_THRESHOLD and not cached:
    # if True:
    if num_outliers <= MAX_OUTLIER_THRESHOLD:
        # generate the PDF
        pdf = run_cmd(["dot", "-Tpdf", dot_file], verbose=False)

        # save the pdf file
        with open(pdf_file, 'w') as fout:
            fout.write(pdf)

    # print some summary stats
    print "{padding}{tree} nodes={nodes}\t admix={admix}\t outliers={out}\t worst={worst}\t {name}".format(
        padding="  "*depth, name=graph_name, tree=newick.ljust(80), nodes=num_nodes, admix=num_admix,
        out=len(outliers), worst=worst_fstat[-1])

    return new_tree, outliers


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
        children = [(child_node.tag, export_newick_tree(child_node)) for child_node in parent_node]
        children.sort()
        tag_name = '' if re.match('n[0-9]+|R', parent_node.tag) else parent_node.tag
        return '(' + ','.join(node for tag, node in children) + ')%s' % tag_name


def run_analysis(all_nodes):
    """
    Build and test all possible trees and graphs
    """

    # setup a simple 2-node tree
    root_node = ElemTree.Element(ROOT)
    root_tree = ElemTree.ElementTree(root_node)

    ElemTree.SubElement(root_node, OUT)
    ElemTree.SubElement(root_node, all_nodes.pop(0))

    # recursively add all the other nodes
    recurse_tree(root_tree, all_nodes[0], all_nodes[1:])

# perform the analysis
run_analysis(NODES)

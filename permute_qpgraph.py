#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Build all possible graphs using a Randomised Stepwise Addition Order Algorithm w/ Branch and Bound.
# Usage...
# python -u permute_qpgraph.py 1> permute-std.log 2> permute-err.log

import xml.etree.ElementTree as ElemTree
import re
import sys
import warnings
import csv

import numpy as np

# import the clustering libraries
from scipy.cluster.hierarchy import linkage, fcluster

# use the Pathos library for improved multi-processing
import pathos.multiprocessing as mp

# import the custom modules
from pipeline_utils import *

from itertools import izip
from cStringIO import StringIO
from Bio import Phylo

with warnings.catch_warnings():
    # dirty hack to suppress warnings from graph_tool
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    from graph_tool.all import *

# TODO improve parsimony...
# on first pass, only allow non-admix insertion
# if can't be added, then send to back of list
# if it fails admix insertion (when it's turn comes up again)
# then throw


class PermuteQpgraph:

    # shoud we use multi-threading to speed up the graph search
    MULTITHREADED_SEARCH = True

    # how many outliers should we allow before pruning a branch in graph space
    MAX_OUTLIER_THRESHOLD = 0

    # print PDFs for graphs with (N - offset) nodes
    REMAINING_PRINT_OFFSET = 0

    def __init__(self, par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive, verbose):
        """
        Initialise the object attributes 
        """
        self.par_file = par_file
        self.dot_path = dot_path
        self.pdf_path = pdf_path
        self.verbose = verbose

        # should we try all possible graphs, or should we stop when we find something reasonable
        self.exhaustive_search = exhaustive

        # open the file for writing
        self.log_handle = open(log_file, 'a')

        if outgroup in nodes:
            nodes.remove(outgroup)

        self.nodes = nodes
        self.outgroup = outgroup

        self.root_node = 'R'
        self.problem_nodes = []
        self.tested_graphs = set()
        self.solutions = set()

    def log(self, message):
        """
        Handle message logging to file/stdout. 
        """
        # send message to the log file
        print >> self.log_handle, message
        self.log_handle.flush()

        if self.verbose:
            # echo to stdout
            print message
            sys.stdout.flush()

    def recurse_tree(self, root_tree, new_tag, remaining, depth=0):
        """
        Permute all possible new trees, by adding the new node to all branches.
    
        If no resulting tree passes the outlier threshold then try adding the node to all possible pairs of branches.
        """
        new_trees = []

        # get all the nodes in the tree (skip the outgroup)
        target_nodes = [node for node in root_tree.findall('.//*') if node.tag != self.outgroup]

        # add the new node to every branch in the tree
        for target_node in target_nodes:

            # clone the current tree and add the new node
            new_tree = copy.deepcopy(root_tree)
            self.insert_node(new_tree, target_node, new_tag)
            new_trees.append(new_tree)

        # test all the trees
        results = self.test_trees(new_trees, depth)

        # process the results
        node_placed = self.check_results(results, remaining, depth)

        # test all the admixture possibilities
        if not node_placed:

            admix_trees = []

            # permute all the two parent admixture possibilities
            pairs = list(itertools.combinations(target_nodes, 2))

            for target1, target2 in pairs:
                # skip duplicate targets (this happens when there is already an admixture node in the tree)
                if target1.tag == target2.tag:
                    continue

                # clone the current tree
                new_tree = copy.deepcopy(root_tree)

                # make a new intermediate node
                admix_label = self.new_label(new_tree, admix=True)

                # add two admix nodes as the children of both targets
                admix_nodes = [
                    self.insert_node(new_tree, target1, admix_label, attrs={'internal': '1', 'admix': '1', 'side': 'l'}),
                    self.insert_node(new_tree, target2, admix_label, attrs={'internal': '1', 'admix': '1', 'side': 'r'})
                ]

                # choose the actual parent based on the sort order of the tag name (needed for unique tree hashing)
                admix_node = admix_nodes[0] if target1.tag < target2.tag else admix_nodes[1]

                # add the new node as the child of the preferred admix node
                self.insert_node(new_tree, admix_node, new_tag, append=True)

                admix_trees.append(new_tree)

            # test all the admixture trees
            results = self.test_trees(admix_trees, depth)

            # process the results
            node_placed = self.check_results(results, remaining, depth)

        if not node_placed:

            # we could not place the node via either method :(
            if new_tag not in self.problem_nodes and remaining and not self.exhaustive_search:
                self.log("WARNING: Unable to place node '%s' at this time." % new_tag)

                self.problem_nodes.append(new_tag)

                # add the problem node to end of the list, as we may be able to add it later on
                remaining.append(new_tag)

                # try and add the other nodes
                self.recurse_tree(root_tree, remaining[0], remaining[1:], depth)

            else:
                raise NodeUnplaceable("ERROR: Cannot place node '%s' in the graph." % new_tag)

    def test_trees(self, new_trees, depth):
        """
        Run qpGraph on a list of trees
        """
        if self.MULTITHREADED_SEARCH:
            # we need to buffer the results to use multi-threading
            pool = mp.ProcessingPool(MAX_CPU_CORES)
            results = pool.map(self.run_qpgraph, itertools.izip(new_trees, itertools.repeat(depth)))
        else:
            # test the trees without multi-threading
            results = []
            for new_tree in new_trees:
                result = self.run_qpgraph((new_tree, depth))
                results.append(result)

        return results

    def check_results(self, results, remaining, depth):
        """
        Check the results from qpGraph
        """
        # were we able to place the new node
        placed_node = False

        for new_tree, outliers, graph_name in results:

            # add this graph to the list of those we've tested
            self.tested_graphs.add(graph_name)

            # did our new trees pass the threshold
            if len(outliers) <= self.MAX_OUTLIER_THRESHOLD:

                # recursively add any remaining nodes
                if remaining:
                    self.recurse_tree(new_tree, remaining[0], remaining[1:], depth + 1)
                else:
                    self.log("SUCCESS: Placed all nodes on a graph without outliers!")

                    # add this graph to the list of solutions
                    self.solutions.add(graph_name)

                # we successfully placed the new node!
                placed_node = True

        return placed_node

    def insert_node(self, new_tree, target_node, new_tag, attrs=None, append=False):
        """
        Helper function to add a new node on the branch leading to the target node.
        """
        # get the target node in the new tree
        target_xpath = './/' + target_node.tag

        if target_node.get('side'):
            target_xpath += '[@side="%s"]' % target_node.get('side')

        target_node = new_tree.find(target_xpath)

        if append:
            # append the new node directly to the target
            parent_node = target_node

        else:
            # get the parent of the target
            parent_node = new_tree.find(target_xpath + '/..')

            # does the target node have a sibling
            if len(parent_node) > 1:
                label = self.new_label(new_tree)
                parent_node.remove(target_node)

                # add an intermediate node, to act as the parent for the new node
                parent_node = ElemTree.SubElement(parent_node, label)
                parent_node.set('internal', '1')

                # re add the target node
                # ElemTree.SubElement(parent_node, target_node.tag)
                parent_node.append(target_node)

        # add the new node as a sibling to the target
        new_node = ElemTree.SubElement(parent_node, new_tag)

        if attrs:
            # add any node attributes
            for key, value in attrs.iteritems():
                new_node.set(key, value)

        return new_node

    def run_qpgraph(self, args):
        """
        Run qpGraph on the given tree 
        """
        # extract the tuple of arguments
        new_tree, depth = args

        # convert the tree to newick format
        newick = self.print_newick_tree(new_tree)

        # get unique names for the output files
        graph_name = self.hash_text(newick)
        grp_file = self.dot_path + '-{name}.graph'.format(name=graph_name)
        dot_file = self.dot_path + '-{name}.dot'.format(name=graph_name)
        log_file = self.dot_path + '-{name}.log'.format(name=graph_name)
        xml_file = self.dot_path + '-{name}.xml'.format(name=graph_name)

        try:
            # if the log file exists then we've run the analysis already
            with open(log_file, 'r') as fin:
                log = fin.read()

        except IOError:
            # save the xml file
            new_tree.write(xml_file)

            # convert the tree to qpGraph format
            graph = self.export_qpgraph(new_tree)

            # save the graph file
            with open(grp_file, 'w') as fout:
                fout.write(graph)

            # run qpGraph
            log = run_cmd(["qpGraph", "-p", self.par_file, "-g", grp_file, "-d", dot_file], verbose=False)

            # save the log file
            with open(log_file, 'w') as fout:
                fout.write(log)

        # parse the log and extract the outliers
        outliers, worst_fstat = self.extract_outliers(log.splitlines())

        # count the leaf nodes
        all_nodes = new_tree.findall('.//*')
        num_nodes = len([node for node in all_nodes if node.get('internal') != '1'])
        num_admix = len([node for node in all_nodes if node.get('admix') == '1']) / 2
        num_outliers = len(outliers)

        # only print PDFs for graphs that pass the threshold
        if num_outliers <= self.MAX_OUTLIER_THRESHOLD and num_nodes > (len(self.nodes) - self.REMAINING_PRINT_OFFSET):

            # embed some useful metadata info in the PDF name
            pdf_file = self.pdf_path + '-n{nodes}-o{out}-a{admix}-{name}.pdf'.format(nodes=num_nodes,
                                                                                     out=num_outliers,
                                                                                     admix=num_admix,
                                                                                     name=graph_name)

            # pretty print the qpGraph dot file
            pprint_qpgraph(dot_file, pdf_file)

        # output some summary stats
        self.log("{padding}{tree} \tnodes={nodes}\t admix={admix}\t outliers={out}\t worst={worst}\t {name}".format(
            padding="  "*depth, name=graph_name, tree=newick.ljust(80), nodes=num_nodes, admix=num_admix,
            out=len(outliers), worst=worst_fstat[-1]))

        return new_tree, outliers, graph_name

    @staticmethod
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

    def export_qpgraph(self, root_tree):
        """
        Convert the ElementTree into qpGraph format
        """
        # clone the tree because this process is destructive
        local_tree = copy.deepcopy(root_tree)

        graph = "root\t{root}\n".format(root=self.root_node)

        # get all the nodes
        nodes = local_tree.findall('.//*')

        for node in nodes:
            # skip non-leaf nodes
            if len(node) == 0 and node.get('admix') != '1':
                graph += "label\t{node}\t{node}\n".format(node=node.tag)

        # build the list of edges
        graph += self.export_qpgraph_node(local_tree)

        return graph

    def export_qpgraph_node(self, root_tree, parent_node=None):
        """
        Recursively export all the edges in the graph
        """
        graph = ""

        if parent_node is None:
            parent_node = root_tree.getroot()

        for child_node in list(parent_node):

            if child_node.get('printed') != '1':

                # is this an admixture node or a normal node
                matches = root_tree.findall('.//' + child_node.tag + '/..')

                if len(matches) > 1:
                    # admixture branch
                    parent1, parent2 = matches

                    middle1 = child_node.tag + 'a'
                    middle2 = child_node.tag + 'b'
                    code1 = self.hash_text(middle1)
                    code2 = self.hash_text(middle2)

                    # don't admix from a bifurcating node; intermediate nodes to accommodate drift
                    graph += "edge\t{code}\t{parent}\t{middle}\n".format(code=code1, parent=parent1.tag, middle=middle1)
                    graph += "edge\t{code}\t{parent}\t{middle}\n".format(code=code2, parent=parent2.tag, middle=middle2)

                    # now admix from the two middle nodes
                    graph += "admix\t{child}\t{parent1}\t{parent2}\t50\t50\n".format(parent1=middle1, parent2=middle2,
                                                                                     child=child_node.tag)

                    # flag both nodes so we don't export them twice
                    parent1.find(child_node.tag).set('printed', '1')
                    parent2.find(child_node.tag).set('printed', '1')

                else:
                    # regular branch
                    code = self.hash_text(child_node.tag)
                    graph += "edge\t{code}\t{parent}\t{child}\n".format(code=code, parent=parent_node.tag,
                                                                        child=child_node.tag)

            # leaf nodes
            if len(child_node) > 0:
                # now convert the children
                graph += self.export_qpgraph_node(root_tree, child_node)

        return graph

    @staticmethod
    def hash_text(text, length=7):
        """
        Generate a unique key by hashing a string
        """
        return hashlib.sha1(text).hexdigest()[0:length]

    @staticmethod
    def new_label(root_tree, admix=False):
        """
        Return a new label for a node
        """
        all_nodes = root_tree.findall('.//*')

        if admix:
            num = len([node for node in all_nodes if node.get('admix') == '1']) / 2
        else:
            num = len([node for node in all_nodes if node.get('internal') == '1' and node.get('admix') != '1'])

        return '{pre}{num}'.format(pre='a' if admix else 'n', num=num+1)

    def print_newick_tree(self, root_tee):
        """
        Convert an ElementTree into a ladderized Newick tree.
        """
        newick = self.export_newick_tree(root_tee.getroot())

        # load into Phylo so we can sort the tree (i.e. ladderize)
        tree = Phylo.read(StringIO(newick), 'newick')
        tree.ladderize()

        # export the tree back to a string
        fout = StringIO()
        Phylo.write(tree, fout, 'newick')

        return fout.getvalue().replace(':0.00000', '').strip()

    def export_newick_tree(self, parent_node):
        """
        Convert an ElementTree tree into Newick format
        """
        if len(parent_node) == 0:
            return parent_node.tag
        else:
            children = [(child_node.tag, self.export_newick_tree(child_node)) for child_node in parent_node]
            children.sort()
            tag_name = '' if re.match('n[0-9]+|R', parent_node.tag) else parent_node.tag
            return '(' + ','.join(node for tag, node in children) + ')%s' % tag_name

    def find_graph(self):
        """
        Build and test all possible trees and graphs
        """

        self.log('INFO: Starting list %s' % self.nodes)

        # setup a simple 2-node tree
        root_node = ElemTree.Element(self.root_node)
        root_tree = ElemTree.ElementTree(root_node)

        ElemTree.SubElement(root_node, self.outgroup)
        ElemTree.SubElement(root_node, self.nodes[0])

        # recursively add all the other nodes
        self.recurse_tree(root_tree, self.nodes[1], self.nodes[2:])


def permute_qpgraph(par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive=False, verbose=False):
    """
    Find the best fitting graph for a given set of nodes, by permuting all possible graphs.
    """

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    qp = PermuteQpgraph(par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive, verbose)

    # get all the permutations of possible node orders
    all_nodes_perms = list(itertools.permutations(nodes, len(nodes)))

    qp.log("INFO: There are %s possible starting orders for the given nodes." % len(all_nodes_perms))
    qp.log("INFO: Performing %s search." % ("an exhaustive" if qp.exhaustive_search else "a heuristic"))

    # keep looping until we find a solution, or until we've exhausted all possible starting orders
    while not qp.solutions or qp.exhaustive_search:

        try:
            # find the best fitting graph for this starting order
            qp.find_graph()

        except NodeUnplaceable as error:
            # log the error
            qp.log(error)

        try:
            # try starting with a different node order
            qp.nodes = list(all_nodes_perms.pop())

        except IndexError:
            # we've run out of node orders to try
            if not qp.solutions:
                qp.log("ERROR: Cannot resolve the graph from any permutation of the given nodes.")

            break

    qp.log("FINISHED: Found %s unique solution(s) from a total of %s unique graphs!" %
           (len(qp.solutions), len(qp.tested_graphs)))

    return len(qp.solutions) > 0


class NodeUnplaceable(Exception):
    """
    Node cannot be placed in the graph without exceeding outlier threshold
    """
    pass


def parse_dot_file(path):
    """
    The graph-tool library doesn't like the header attributes used by qpGraph, so we need to filter them out
    """
    with open(path, 'r') as fin:
        rows = fin.readlines()

    # exclude lines 2-4, which contain the problematic metadata
    text = "".join(rows[:1] + rows[5:])

    return StringIO(text)


def find_clusters(graph_names, pdf_file, csv_file):

    from scipy.cluster.hierarchy import linkage

    graphs = []

    # instantiate all the graph objects
    for graph_name in graph_names:
        # TOOD self.dot_path
        dot_file = dot_path + '-{name}.dot'.format(name=graph_name)
        graph = load_graph(parse_dot_file(dot_file), fmt='dot')
        graphs.append(graph)

    # how many graphs are we comparing
    size = len(graph_names)

    # initialise a distance matrix
    dist_matrix = np.zeros([size, size])

    # populate the distance matrix
    for i in range(1, size):
        for j in range(i):
            # calculate the distance scores between graph pairs (scores are not symmetric; i.e. A->B != B->A)
            d1 = 1 - similarity(graphs[i], graphs[j])
            d2 = 1 - similarity(graphs[j], graphs[i])

            # enforce symmetry by taking the max distance
            dist_matrix[i, j] = dist_matrix[j, i] = max(d1, d2)

    # calculate the hierarchical clusters, using Ward's minimum variance method
    # see https://en.wikipedia.org/wiki/Ward%27s_method
    linkage = linkage(dist_matrix, method='ward')

    # print a dendrogram of the clusters
    # pprint_dendrogram(
    #     linkage,
    #     truncate_mode='lastp',
    #     p=10,
    #     leaf_rotation=90.,
    #     leaf_font_size=12.,
    #     show_contracted=True,
    #     pdf=pdf_file,
    #     # max_d=20,  # plot a horizontal cut-off line
    # )

    # automatically assign graphs to clusters
    # https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/#Inconsistency-Method
    clusters = fcluster(linkage, 8, criterion='inconsistent', depth=10)

    print "Found %s clusters" % len(set(clusters))

    with open(csv_file, 'wb') as fout:
        csv_writer = csv.writer(fout)
        csv_writer.writerow(['Graph', 'Cluster'])
        for graph, cluster in izip(graph_names, clusters):
            csv_writer.writerow([graph, cluster])

    return clusters


if __name__ == "__main__":

    # simulated test data...
    # par_file = 'permute/simulated.par'
    # log_file = 'permute/simulated.log'
    # dot_path = 'permute/graphs/sim'
    # pdf_path = 'permute/pdf/sim'
    # nodes = ['A', 'B', 'C', 'X']
    # outgroup = 'Out'

    # if len(sys.argv) != 3:
    #     print "Error: required params"
    #     quit()
    #
    # group = sys.argv[1]
    # dataset = sys.argv[2]
    #
    # nodes = GROUPS[dataset][group]
    # outgroup = OUTGROUP_POP[group] if group in OUTGROUP_POP else OUTGROUP_POP[dataset]
    #
    # par_file = 'qpgraph/{0}.{1}.permute.par'.format(group, dataset)
    # log_file = 'qpgraph/{0}.{1}.permute.log'.format(group, dataset)
    # dot_path = 'qpgraph/{0}.permute'.format(dataset)
    # pdf_path = 'pdf/{0}.{1}.qpg-permute'.format(group, dataset)


    # permute_qpgraph(par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive=True, verbose=True)

    # TODO fix me
    import glob
    files = glob.glob('pdf/graph-pops2.merged_v2_TV_laurent.qpg-permute-*')
    graph_names = [re.search(r'a[0-9]-(.+).pdf', file).group(1) for file in files]

    dot_path = 'qpgraph/merged_v2_TV_laurent.permute'

    # graph_names = ['03f8f20', '04415fe', '074f42f', '0abe334', '0bd6ab7', '0cc8dc9', '0eb6d63', '1169e75', '15ff806',
    #                '16a4713', '19ed6e8', '1b58e0b', '1e3e40e', '1e5f242', '1edf0bb', '1fa4501', '204a55b', '20edf7c',
    #                '2280052', '252a051', '256b3f9', '29925e2', '2a879fc', '2e5374a', '2e6f645', '30f25d3', '31790c8',
    #                '31c40c7', '3209cf8', '322bc1b', '3274050', '35b7f88', '3860a94', '3ba1bdf', '3e309e8', '3e6cd23',
    #                '401685b', '405ff67', '410ecac', '43ec82c', '449a55f', '46c0f99', '4add777', '4c093f0', '4c8ccfa',
    #                '5055248', '51eb595', '52443ef', '5342b7f', '5364738', '537fb37', '53ceace', '56a7d9a', '5ba7166',
    #                '5c3c62a', '5df249f', '649e55c', '64ebc31', '68029e4', '6946338', '6b7545a', '6bd8707', '6d75411',
    #                '6d95fca', '725683f', '725a0f0', '72b9c75', '769a660', '76da55b', '779ab27', '77de3f5', '78df87f',
    #                '7d3fbad', '7e6156c', '7f0646c', '801c201', '8144745', '82ea5f2', '83a8650', '863a17b', '86470cc',
    #                '87369c9', '8a40fb1', '8b0ebca', '9773ea2', '98617f0', '995962b', '9ae068c', '9eada70', 'a19cbf2',
    #                'a3013ca', 'a3177ad', 'a58b8cd', 'a63dfbf', 'a78ade8', 'a85f651', 'aa24aee', 'abe0beb', 'ac481d3',
    #                'ae3d738', 'af2cc15', 'b2b46e0', 'b2d5e17', 'b443fe6', 'b564cc5', 'b5dfeb9', 'b767d73', 'b82be21',
    #                'b85acd2', 'b86960d', 'b988ce4', 'ba29db9', 'bcbb494', 'c0f9e78', 'c35b8da', 'c3c7908', 'c82d89f',
    #                'c836ff0', 'c8f0147', 'ca128bd', 'ca82eef', 'cea743e', 'd1edd48', 'd35ace2', 'd937116', 'd97a1b5',
    #                'da42c84', 'dae410e', 'db585ea', 'dbb73f9', 'dbf0daa', 'dd3e2a1', 'ddd2b7d', 'de54ee0', 'e35cad4',
    #                'e5c3cc6', 'eb9a011', 'ec3e11f', 'ef3d704', 'ef9c447', 'f019c0f', 'f0468a5', 'f57d5b5', 'f98e287',
    #                'fdbe395', 'ffbcae6']

    find_clusters(graph_names, pdf_file='clusters.pdf', csv_file='clusters.csv')

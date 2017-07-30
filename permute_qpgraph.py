#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Build all possible graphs using a Randomised Stepwise Addition Order Algorithm w/ Branch and Bound.
# Usage...
# python -u permute_qpgraph.py 1> permute-std.log 2> permute-err.log

import xml.etree.ElementTree as ElemTree
import re
import sys
import csv
import glob

# TODO this throws a warning...
# import matplotlib
# matplotlib.use('Agg')

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

from graph_tool import *
from graph_tool.topology import *

# TODO improve parsimony...
# on first pass, only allow non-admix insertion
# if can't be added, then send to back of list
# if it fails admix insertion (when it's turn comes up again)
# then throw


class PermuteQpgraph:

    # should we use multi-threading to speed up the graph search
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

        If that still doesn't work, then try inserting the node as a child of an existing branch, but also as an admix
        parent of an existing node. For example, (A,(B,X)) + C -> ((A,(B,(X)a1)),(C)a1)
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

        # test all the admixture child possibilities
        if not node_placed:

            admix_child_trees = []

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

                admix_child_trees.append(new_tree)

            # test all the admixture trees
            results = self.test_trees(admix_child_trees, depth)

            # process the results
            node_placed = self.check_results(results, remaining, depth)

        # test all the admixture parent possibilities
        if not node_placed:

            admix_parent_trees = []

            # none of the normal insertions worked, so lets try them all again but this time lets allow an existing
            # node to be removed and reinserted as the admixed child of its current parent and the new node
            for outlier_tree in new_trees:

                # get the new node and its parent
                new_node = outlier_tree.find('.//' + new_tag)
                parent_node = outlier_tree.find('.//' + new_tag + '/..')

                # try all possible new admixed children
                for target_node in target_nodes:

                    # construct an xpath query to find the target in the current tree
                    target_xpath = './/' + target_node.tag

                    if target_node.get('side'):
                        target_xpath += '[@side="%s"]' % target_node.get('side')

                    # skip targets which are ancestral to the new node
                    if outlier_tree.find(target_xpath + '//' + new_node.tag) is not None:
                        continue

                    # clone the current tree
                    new_tree = copy.deepcopy(outlier_tree)

                    # get the target and it's parent in this copy of the tree
                    target_node = new_tree.find(target_xpath)
                    target_parent = new_tree.find(target_xpath + '/..')

                    # make a new intermediate node
                    admix_label = self.new_label(new_tree, admix=True)

                    # remove the target from the tree
                    target_parent.remove(target_node)

                    # add two admix nodes as the children of both targets
                    admix_nodes = [
                        self.insert_node(new_tree, parent_node, admix_label, attrs={'internal': '1', 'admix': '1', 'side': 'l'}),
                        self.insert_node(new_tree, target_parent, admix_label, attrs={'internal': '1', 'admix': '1', 'side': 'r'})
                    ]

                    # choose the actual parent based on the sort order of the tag name (needed for unique tree hashing)
                    admix_node = admix_nodes[0] if parent_node.tag < target_parent.tag else admix_nodes[1]

                    # re add the target node
                    admix_node.append(target_node)

                    # add the new tree to the list of trees to test
                    admix_parent_trees.append(new_tree)

            # test all the admixture parent trees
            results = self.test_trees(admix_parent_trees, depth)

            # process the results
            node_placed = self.check_results(results, remaining, depth)

        if not node_placed:

            self.log("WARNING: Unable to place node '%s'." % new_tag)
            raise NodeUnplaceable("ERROR: Cannot place node '%s' in the graph." % new_tag)

            # # we could not place the node via either method :(
            # if new_tag not in self.problem_nodes and remaining and not self.exhaustive_search:
            #     self.log("WARNING: Unable to place node '%s' at this time." % new_tag)
            #
            #     self.problem_nodes.append(new_tag)
            #
            #     # add the problem node to end of the list, as we may be able to add it later on
            #     remaining.append(new_tag)
            #
            #     # try and add the other nodes
            #     self.recurse_tree(root_tree, remaining[0], remaining[1:], depth)
            #
            # else:
            #     raise NodeUnplaceable("ERROR: Cannot place node '%s' in the graph." % new_tag)

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
        newick = fout.getvalue()

        # remove the branch lenghs
        newick = newick.replace(':0.00000', '').strip()

        # get the order of admix nodes in the tree
        order = list(OrderedDict.fromkeys(re.findall('a\d+', newick)))

        # normalise the node numbering
        for i, old in enumerate(order):
            newick = newick.replace(old, 'n%s' % (i+1))

        # replace n0 with a0 (to preseve the existing chache)
        newick = re.sub(r'n([0-9]+)', r'a\1', newick)

        return newick

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


class NodeUnplaceable(Exception):
    """
    Node cannot be placed in the graph without exceeding outlier threshold
    """
    pass


def permute_qpgraph(par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive=False, verbose=False):
    """
    Find the best fitting graph for a given set of nodes, by permuting all possible graphs.
    """

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    pq = PermuteQpgraph(par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive, verbose)

    pq.log("INFO: Performing %s search." % ("an exhaustive" if pq.exhaustive_search else "a heuristic"))

    try:
        # find the best fitting graph
        pq.find_graph()

    except NodeUnplaceable as error:
        # log the error
        pq.log(error)

        # we didn't find a solution
        if not pq.solutions:
            pq.log("ERROR: Cannot resolve the graph all the given nodes.")

    pq.log("FINISHED: Found %s unique solution(s) from a total of %s unique graphs!" %
           (len(pq.solutions), len(pq.tested_graphs)))

    return pq.solutions


class ClusterQpgraph():

    # should we use multi-threading to speed up the graph search
    MULTITHREADED_SEARCH = True

    def __init__(self, graph_names, log_file, dot_path, csv_file, mtx_file, verbose):
        """
        Initialise the object attributes
        """
        self.graph_names = graph_names
        self.log_file = log_file
        self.dot_path = dot_path
        self.csv_file = csv_file
        self.mtx_file = mtx_file
        self.verbose = verbose

        self.graphs = []

        self.verbose = True

        # open the file for writing
        self.log_handle = open(log_file, 'a')

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

    @staticmethod
    def parse_dot_file(path):
        """
        The graph-tool library doesn't like the header attributes used by qpGraph, so we need to filter them out
        """
        with open(path, 'r') as fin:
            rows = fin.readlines()

        # exclude lines 2-4, which contain the problematic metadata
        text = "".join(rows[:1] + rows[5:])

        return StringIO(text)

    def calculate_distance(self, args):
        """
        Calculate the similarity distance for two graphs.

        See https://graph-tool.skewed.de/static/doc/topology.html#graph_tool.topology.similarity
        """

        # extract the tuple of arguments
        i, j = args

        # TODO check how node names effect this!
        # calculate the distance scores between graph pairs (scores are not symmetric; i.e. A->B != B->A)
        d1 = similarity(self.graphs[i], self.graphs[j], distance=True)
        d2 = similarity(self.graphs[j], self.graphs[i], distance=True)

        # enforce symmetry in the matrix by taking the max distance
        dist = max(d1, d2)

        return i, j, dist

    def build_matrix(self):
        """
        Build a symmetrical distance matrix for all graphs.
        """

        # instantiate all the graph objects
        for graph_name in self.graph_names:
            dot_file = self.dot_path + '-{name}.dot'.format(name=graph_name)
            graph = load_graph(self.parse_dot_file(dot_file), fmt='dot')
            self.graphs.append(graph)

        # how many graphs are we comparing
        size = len(self.graph_names)

        # initialise the distance matrix
        dist_matrix = np.zeros([size, size])

        # get all the i,j pairs for one diagonal half
        idxs = [(i, j) for i in range(1, size) for j in range(i)]

        self.log("INFO: Calculating distance matrix for %s graph pairs" % len(idxs))

        if self.MULTITHREADED_SEARCH:
            # we need to buffer the results to use multi-threading
            pool = mp.ProcessingPool(MAX_CPU_CORES)
            results = pool.map(self.calculate_distance, idxs)
        else:
            # compute distances without multi-threading
            results = []
            for i, j in idxs:
                result = self.calculate_distance((i, j))
                results.append(result)

        # populate the distance matrix
        for i, j, dist in results:
            dist_matrix[i, j] = dist_matrix[j, i] = dist

        # save the matrix
        np.save(self.mtx_file, dist_matrix)

        return dist_matrix

    def get_matrix(self):
        """
        Load the distance matix from file, or build it if necessary.
        """
        try:
            # load the distance matrix from file
            dist_matrix = np.load(self.mtx_file)

            self.log("INFO: Loaded distance matrix from file %s" % self.mtx_file)

        except IOError:
            # file doesn't exist, so build it
            dist_matrix = self.build_matrix()

        return dist_matrix


def cluster_qpgraph(graph_names, dot_path, log_file, pdf_file, csv_file, mtx_file, verbose=False):
    """
    Compare all fitting graphs and compute the number of clusters.
    """

    # clean up the log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # instantiate the class
    cq = ClusterQpgraph(graph_names, log_file, dot_path, csv_file, mtx_file, verbose)

    cq.log("INFO: There are %s graphs to compare" % len(set(graph_names)))

    # get the distance matrix
    dist_matrix = cq.get_matrix()

    cq.log("INFO: Calculating the hierarchical clusters (linkage matrix)")

    # calculate the hierarchical clusters, using Ward's minimum variance method
    # https://en.wikipedia.org/wiki/Ward%27s_method
    Z = linkage(dist_matrix, method='ward')

    # print a dendrogram of the clusters
    pprint_dendrogram(Z, truncate_mode='lastp', p=10, leaf_rotation=90.,
                      leaf_font_size=12., show_contracted=True, pdf=pdf_file)

    cq.log("INFO: Printed hierarchical clustering dendrogram %s" % pdf_file)

    # automatically assign graphs to clusters
    # https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/#Inconsistency-Method
    clusters = fcluster(Z, t=10, criterion='inconsistent', depth=10)

    cq.log("INFO: Found %s clusters using inconsistency criterion (t=%s)" % (len(set(clusters)), 10))

    with open(csv_file, 'wb') as fout:
        csv_writer = csv.writer(fout)
        csv_writer.writerow(['Graph', 'Cluster'])
        for graph, cluster in izip(graph_names, clusters):
            csv_writer.writerow([graph, cluster])

    cq.log("INFO: Saved clusters to file %s" % csv_file)


if __name__ == "__main__":

    from time import time
    from datetime import timedelta

    start = time()

    # simulated test data...
    nodes = ['A', 'X', 'B', 'C']
    outgroup = 'Out'
    par_file = 'permute/simulated.par'
    log_file = 'permute/simulated.log'
    dot_path = 'permute/graphs/sim'
    pdf_path = 'permute/pdf/sim'

    # if len(sys.argv) != 3:
    #     print "Error: required params"
    #     quit()
    #
    # group = sys.argv[1]
    # dataset = sys.argv[2]
    #
    # dot_path = 'qpgraph/dot/{0}.permute'.format(dataset)
    #
    # # ---------------------
    # # -- PERMUTE_QPGRAPH --
    # # ---------------------
    # nodes = GROUPS[dataset][group]
    # outgroup = OUTGROUP_POP[group] if group in OUTGROUP_POP else OUTGROUP_POP[dataset]
    # par_file = 'qpgraph/{0}.{1}.permute.par'.format(group, dataset)
    # log_file = 'qpgraph/{0}.{1}.permute.log'.format(group, dataset)
    # pdf_path = 'pdf/{0}.{1}.qpg-permute'.format(group, dataset)

    permute_qpgraph(par_file, log_file, dot_path, pdf_path, nodes, outgroup, exhaustive=True, verbose=True)

    # ---------------------
    # -- CLUSTER_QPGRAPH --
    # ---------------------
    # log_file = 'qpgraph/{0}.{1}.cluster.log'.format(group, dataset)
    # csv_file = 'qpgraph/{0}.{1}.cluster.csv'.format(group, dataset)
    # mtx_file = 'qpgraph/{0}.{1}.cluster.npy'.format(group, dataset)
    # pdf_file = 'pdf/{0}.{1}.qpg-cluster.pdf'.format(group, dataset)
    #
    # # find all the PDFs, and extract the graph names
    # files = glob.glob('pdf/{0}.{1}.qpg-permute-*'.format(group, dataset))
    # graph_names = [re.search(r'a[0-9]-(.+).pdf', file).group(1) for file in files]
    #
    # cluster_qpgraph(graph_names, dot_path, log_file, pdf_file, csv_file, mtx_file, verbose=True)

    print "INFO: Execution took: %s" % timedelta(seconds=time()-start)
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyparsing, re, gzip, sys, json, hashlib

if len(sys.argv) != 2:
    print "Error: Missing treeout.gz parameter"
    quit()

tree_path = sys.argv[1]

labels = {}


def new_label():
    """
    Return a new label for a node
    """
    new_label.n += 1
    return 'n{}'.format(new_label.n)
new_label.n = 0


def treemix_tree(nested):
    """
    Recursively resolve the NODE:DRIFT notation used by treemix
    """
    nodes = []
    children = prev_item = None

    for item in nested:
        if isinstance(item, list):
            children = treemix_tree(item)
        else:
            m = re.match('^([A-Za-z0-9]+)?:([0-9.e-]+)$', item)
            if m.group(1):
                # leaf node (ie. LABEL:DRIFT)
                hash = hash_node(item)
                labels[hash] = m.group(1)
                node = (hash, m.group(2), [])
            else:
                # internal node (i.e. :DRIFT)
                hash = hash_node([prev_item, item])
                labels[hash] = new_label()
                node = (hash, m.group(2), children)
            nodes.append(node)
        prev_item = item
    return nodes


def print_treemix_tree(tree, parent, migs):
    """
    Output the treemix tree as a series of edges and admix nodes
    """
    for node in tree:
        hash, drift, children = node

        if hash not in migs:
            print "edge\t%s\t%s\t%s" % (hash, parent, labels[hash])
        else:
            migrant, percent = migs[hash]
            print "admix\t%s\t%s\t%s\t%s\t%s" % (labels[hash], parent, labels[migrant], (1-percent), percent)

        if len(children):
            print_treemix_tree(children, labels[hash], migs)


def hash_node(node):
    """
    Hash a treemix node (so we can uniquely identify node positions)
    """
    text = json.dumps(node)
    for item, rep in [('[','('), (']',')'), ('"',''), (' ',''), (',:',':')]:
        text = text.replace(item, rep)
    return hash_text(text if isinstance(node, str) else text[1:-1])


def hash_text(text, len=7):
    """
    Hash a string
    """
    return hashlib.sha1(text).hexdigest()[0:len]


with gzip.open(tree_path, 'r') as fin:

    # fetch the tree topology
    topology = fin.readline()

    # get all the leaf nodes
    leaves = re.findall("([A-Za-z]+):", topology)

    # get a simple nested expression parser
    parser = pyparsing.nestedExpr()

    # use whitespace delimiters
    topology = topology.replace(',', ' ')

    # parse the topology into a nested list
    nested = parser.parseString(topology).asList().pop()

    # parse the list into a tree
    tree = treemix_tree(nested)

    # make a dict of migration branches
    migs = {}
    for line in fin:
        line = line.split()
        percent, parent, child = [line[i] for i in [0, 4, 5]]
        migs[hash_text(child)] = hash_text(parent), float(percent)

    print "root\tR"
    for leaf in leaves:
        print "label\t{label}\t{label}".format(label=leaf)

    print_treemix_tree(tree, 'R', migs)

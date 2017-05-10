#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from Bio import Phylo
import xml.etree.ElementTree as ET
import copy, hashlib

root = 'R'
out = 'O'

nodes = ['A', 'B', 'C', 'X']


def add_node(root_tree, new_node):

    # get all the leaf nodes
    nodes = root_tree.findall('.//*')

    for node in nodes:
        # skip the outgroup and non-leaf nodes
        if node.tag == out or len(node) > 0:
            continue

        # clone the current tree
        new_tree = copy.deepcopy(root_tree)

        ET.dump(new_tree)

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
        ET.SubElement(root_node, new_node)

        yield new_tree

def build_tree(root_tree, unplaced):
    CNTR = 1

    for new_node in unplaced:

        # adding each new node generates (n-1)*2 new trees, where n is the number of nodes already in the tree
        new_trees = add_node(root_tree, new_node)

        # get the remaining nodes
        remaining = list(unplaced)
        remaining.remove(new_node)

        for new_tree in new_trees:
            ET.dump(new_tree)

            print convert_tree(new_tree)

            CNTR += 1

            if CNTR > 2:
                quit()

        #     # TODO export tree to qpgraph format
        #
        #     # TODO run qpgraph
        #
        #     # TODO parse the result, quantify the outliers
        #     # count of outliers
        #     # worst f-stat
        #
        #     # TODO for each new tree that passes threshold, lets add the remaining nodes
            build_tree(new_tree, remaining)


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
    """
    Hash a string
    """
    return hashlib.sha1(text).hexdigest()[0:len]

def new_label():
    """
    Return a new label for a node
    """
    new_label.n += 1
    return 'n{}'.format(new_label.n)
new_label.n = 0

# lets get started
for node in nodes:

    # setup the simplest 2-node tree
    root_node = ET.Element(root)
    root_tree = ET.ElementTree(root_node)

    ET.SubElement(root_node, out)
    ET.SubElement(root_node, node)

    # TODO debugging
    # print "--------"
    # ET.dump(root_tree)

    # get the unplaced nodes
    unplaced = list(nodes)
    unplaced.remove(node)

    # print "Unplaced... %s" % unplaced

    build_tree(root_tree, unplaced)

    break

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from Bio import Phylo
import xml.etree.ElementTree as ET
import copy, hashlib

root = 'R'
out = 'O'

nodes = ['A', 'B', 'C', 'X']


def add_node(sub_tree, new_node):

    # get all the leaf nodes
    nodes = sub_tree.findall('*')

    for node in nodes:
        # skip the outgroup
        if node.tag == out:
            continue

        # clone the current tree
        new_tree = copy.deepcopy(sub_tree)

        # get the parent node
        parent_node = new_tree.find(node.tag)

        # add new node as child
        ET.SubElement(parent_node, new_node)

        yield new_tree


def build_tree(sub_tree, unplaced):

    for new_node in unplaced:

        # adding each new node generates (n-1)*2 new trees, where n is the number of nodes already in the tree
        new_trees = add_node(sub_tree, new_node)

        # get the remaining nodes
        remaining = list(unplaced)
        remaining.remove(new_node)

        for new_tree in new_trees:
            ET.dump(new_tree)

            print convert_tree(new_tree)


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


def convert_tree(sub_tree):
    new_label.n = 0

    graph = "root\t{root}\n".format(root=root)
    graph += "label\t{out}\t{out}\n".format(out=out)

    for node in nodes:
        graph += "label\t{node}\t{node}\n".format(node=node)

    graph += convert_node(sub_tree, sub_tree.tag)

    return graph


def convert_node(parent_node, parent_tag):
    graph = ""

    for child_node in parent_node:

        child_code = hash_text(child_node.tag)

        # leaf nodes
        if len(child_node) == 0:
            graph += "edge\t{code}\t{parent}\t{inter}\n".format(code=child_code, parent=parent_tag, inter=child_node.tag)

        else:
            # insert an intermediate node between the parent and child
            inter = new_label()
            parent_code = hash_text(inter)

            graph += "edge\t{code}\t{parent}\t{inter}\n".format(code=parent_code, parent=parent_tag, inter=inter)
            graph += "edge\t{code}\t{inter}\t{child}\n".format(code=child_code, inter=inter, child=child_node.tag)

            # now convert the children
            graph += convert_node(child_node, inter)

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

# lets get started
for node in nodes:

    # setup the simplest 2-node tree
    root_node = ET.Element(root)
    # root_tree = ET.ElementTree(root_node)
    ET.SubElement(root_node, out)
    ET.SubElement(root_node, node)

    # TODO debugging
    # print "--------"
    # ET.dump(root_tree)

    # get the unplaced nodes
    unplaced = list(nodes)
    unplaced.remove(node)

    # print "Unplaced... %s" % unplaced

    build_tree(root_node, unplaced)

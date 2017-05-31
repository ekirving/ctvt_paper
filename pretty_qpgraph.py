#!/usr/bin/env python
# -*- coding: utf-8 -*-

import graphviz
import re
import csv

# the name of the graph file
dot_file = "permute/_final/0f3755a.dot"
pdf_file = "permute/_final/0f3755a"  # no extension

csv_file = "pop_names.csv"

# load the file
with open(dot_file, 'rU') as fin:
    file = fin.read()

# extract the body contents from the file
body = re.search("{(.+)}", file, re.DOTALL).group(1).strip().split("\n")

# make a new direct graph
dot = graphviz.Digraph(body=body)

# remove the messy graph label
dot.attr('graph', label='')

# set Node defaults
dot.node_attr['shape'] = 'point'
dot.node_attr['fontname'] = 'arial'
dot.node_attr['fontsize'] = '11'

# set Edge defaults
dot.edge_attr['arrowhead'] = 'vee'
dot.edge_attr['fontcolor'] = '#838b8b'  # grey
dot.edge_attr['fontname'] = 'arial'
dot.edge_attr['fontsize'] = '11'

nodes = []

for line in dot:
    # extract the leaf nodes from the body of the graph
    match = re.search("^ +([a-z]+) +\[", line, re.IGNORECASE)
    if match:
        nodes.append(match.group(1))

# import colors from project
with open(csv_file, 'rU') as fin:
    colors = dict((line[0], line[3]) for line in csv.reader(fin))

# set leaf node attributes
for node in nodes:
    dot.node(node, shape='ellipse',
                   color=colors[node],
                   fontcolor=colors[node])

# render the graph
dot.render(pdf_file)
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import graphviz
import re

# import the custom modules
from pipeline_utils import *

# import csv
# import colors from project
# with open('pop_names.csv', 'rU') as fin:
#     colors = dict((line[0], line[3]) for line in csv.reader(fin))

# the base name of the dot file
basename = "permute/_final/0f3755a"

# extract the body contents from the dot file
with open(basename + '.dot', 'rU') as fin:
    body = re.search("{(.+)}", fin.read(), re.DOTALL).group(1).strip().split("\n")

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

# set leaf node attributes
for node in nodes:
    colour = COLOURS.get(POPULATIONS.get(node), DEFAULT_COLOUR)
    dot.node(node, shape='ellipse', color=colour, fontcolor=colour)

# render the graph
dot.render(basename)

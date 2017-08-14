from graph_tool import *
from graph_tool.topology import *
from cStringIO import StringIO

def parse_dot_file(path):
    """
    The graph-tool library doesn't like the header attributes used by qpGraph, so we need to filter them out
    """
    with open(path, 'r') as fin:
        rows = fin.readlines()

    # exclude lines 2-4, which contain the problematic metadata
    text = "".join(rows[:1] + rows[5:])

    return StringIO(text)

code1 = 'bc7e8ff'    # (C,(B,A))
code2 = 'bc7e8ff-2'  # (C,(A,B))

# things that made no difference...
# - branch labels / i.e. drift params
# - order of nodes/branches in the dot file

# things that do make a difference...
# - node names where the chrs are different, or
# - node names where the number decreases
# - n1 == n2 == n3
# - n2 != n1 != n3

# The adjacency similarity is the sum of equal entries in the adjacency matrix, given a vertex ordering determined by
# the vertex labels. In other words, it counts the number of edges which have the same source and target labels in both
# graphs.

# code2 = '8416649'   # (C,(A,X))
# code3 = '0fdeeb2'  # (C,(X,B))
# code4 = '0c82ed9'  # (X,(A,B))
# code5 = '7ef357c'  # (B,(C,A))

code_a = code1
code_b = code2

g1_path = 'permute/graphs/sim-{}.dot'.format(code_a)
g2_path = 'permute/graphs/sim-{}.dot'.format(code_b)

print g1_path
print g2_path

g1 = load_graph(parse_dot_file(g1_path), fmt='dot')
g2 = load_graph(parse_dot_file(g2_path), fmt='dot')

# print g1.vertex_index
# print g1.vertex_properties['vertex_name']
# print g1.vertex_properties['label']
print g1.vertex_index
for key in g1.edges():
    # print key
    print g1.edge(*key)
# print PropertyMap(g1, 'v')


# print g1.vertex_properties['label'].copy('unsigned long')

g1.vertex_properties["foo"] = g1.new_vertex_property('unsigned long', vals='')
print g1.vertex_properties["foo"].value_type()

print g1.list_properties()



label1 = g1.vertex_index
label2 = g2.vertex_index

# print label1.value_type()

print "similarity(%s, %s) = %s " % (code_a, code_b, similarity(g1, g2, distance=False))
print "similarity(%s, %s) = %s " % (code_a, code_b, similarity(g1, g2, label1=label1, label2=label2, distance=False))
print ""
print "similarity(%s, %s) = %s " % (code_b, code_a, similarity(g2, g1, distance=False))

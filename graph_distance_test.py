from graph_distance import *

graph1 = [ ("a","b"), ("b","c"), ("b","d"), ("d","e"), \
           ("e","f"), ("b","f"), ("b","g"), ("f", "g"),
           ("a","g"), ("a","g"), ("c","d"), ("d", "g"),
           ("d","h"), ("aa","h"), ("aa","c"), ("f", "h"), ]

graph2 = [ ("a","b2"), ("b2","c"), ("b2","d"), ("d","e"), \
           ("e","f"), ("b2","f"), ("b2","g"), ("f", "g"),
           ("a","g"), ("a","g"), ("c","d"), ("d", "g"),
           ("d","h"), ("aa","h"), ("aa","c"), ("f", "h"), ]

graph1 = GraphDistance(graph1)
graph2 = GraphDistance(graph2)

distance, graph = graph1.distance_matching_graphs_paths(graph2, use_min=False, store=None)

print distance

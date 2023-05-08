
import numpy as np
from typing import List, Tuple

def is_max(G, probabilities, true_index):
    if probabilities[true_index] == max(probabilities):
        return 1
    else:
        return 0

def dist_from_max(G, probabilities, true_index):
    min_dist = 100000
    eps = 0.000001
    max_val = max(probabilities)
    graph_diameter = np.amax(G)
    for i in range(len(probabilities)):
        if probabilities[i] == max_val:
            min_dist = min(min_dist, G[true_index, i])
    return float(1 - np.exp(-min_dist/(graph_diameter-min_dist+eps)))

def average_dist_from_max(G, probabilities, true_index):
    eps = 0.000001
    max_val = max(probabilities)
    graph_diameter = np.amax(G)
    dists = 0
    count = 0
    for i in range(len(probabilities)):
        if probabilities[i] == max_val:
            value = 1 - np.exp(-G[true_index,i]/(graph_diameter-G[true_index,i]+eps))
            dists += value * probabilities[i]
            count += probabilities[i]
    return float(dists/count)

def average_dist(G, probabilities, true_index):
    eps = 0.000001
    graph_diameter = np.amax(G)
    dists = 0
    count = 0
    for i in range(len(probabilities)):
        value = 1 - np.exp(-G[true_index, i]/(graph_diameter-G[true_index, i]+eps))
        dists += value * probabilities[i]
        count += probabilities[i]
    return float(dists/count)

# def calculate_spanning_graph(G, probabilities):
#     max_val = max(probabilities)
#     graph = []
#     for i in range(len(probabilities)):
#         for j in range(len(probabilities)):
#             if probabilities[i] == max_val and probabilities[j] == max_val:
#                 graph.append((i, j, G[i][j]))
#     return graph

# def kruskal(graph: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
#     sorted_edges = sorted(graph, key=lambda x: x[2])
#     components = {v: k for k, v in enumerate(set([u for u, _, _ in graph] + [v for _, v, _ in graph]))}
#     min_spanning_tree = []
#     for edge in sorted_edges:
#         u, v, w = edge
#         # If the endpoints of the edge belong to different connected components, add the edge to the minimum spanning tree
#         if components[u] != components[v]:
#             min_spanning_tree.append(edge)
            
#             # Merge the connected components
#             old_component = components[v]
#             new_component = components[u]
#             for vertex, component in components.items():
#                 if component == old_component:
#                     components[vertex] = new_component
    
#     return min_spanning_tree

# def calculate_minimum_spanning_graph(G, probabilities):
#     spanning_graph = calculate_spanning_graph(G, probabilities)

#     # floyd-warshall algorithm
#     for k in range(len(spanning_graph)):
#         for i in range(len(spanning_graph)):
#             for j in range(len(spanning_graph)):
#                 spanning_graph[i][j] = min(spanning_graph[i][j], spanning_graph[i][k] + spanning_graph[k][j])
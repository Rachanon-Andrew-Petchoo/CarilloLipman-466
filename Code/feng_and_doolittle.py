from pairwise_alignment import pairwise_align
from progressive_alignment import progressive_alignment
from utils import get_delta_for_max, get_delta_for_min
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

delta_for_max = None
delta_for_min = None

def feng_doolittle_alignment(sequences, alphabet):
    global delta_for_max, delta_for_min
    if delta_for_max == None:
       delta_for_max = get_delta_for_max()
    if delta_for_min == None:
       delta_for_min = get_delta_for_min()
    n = len(sequences)

    # 1. Create distance matrix
    distance_matrix = np.zeros((n,n))
    for i in range(n):
       for j in range(n):
         score = (-1) * pairwise_align(sequences[i], sequences[j], delta_for_max)[0]
         distance_matrix[i][j] = score
         distance_matrix[j][i] = score

    # print(f'Distance matrix: {distance_matrix}')

    # 2. Create graph & MST
    graph = nx.Graph()

    graph.add_nodes_from([ (i, {"alignment": [sequences[i]]}) for i in range(n) ])
    for i in range(n):
      for j in range(i + 1, n):
          weight = distance_matrix[i][j]
          graph.add_edge(i, j, weight=weight)

    mst = nx.minimum_spanning_tree(graph)

    # Visualize the original graph and the MST
    # pos = nx.spring_layout(graph)  # positions for all nodes
    # plt.figure(figsize=(10, 5))

    # plt.subplot(1, 2, 1) # Original graph
    # nx.draw(graph, pos, with_labels=True)
    # labels = nx.get_edge_attributes(graph, 'weight')
    # nx.draw_networkx_edge_labels(graph, pos, edge_labels=labels)
    # plt.title('Original Graph')

    # plt.subplot(1, 2, 2) # MST
    # nx.draw(mst, pos, with_labels=True, edge_color='red')
    # labels = nx.get_edge_attributes(mst, 'weight')
    # nx.draw_networkx_edge_labels(mst, pos, edge_labels=labels)
    # plt.title('Minimum Spanning Tree')

    # plt.show()

    # 3. Align using the guide tree

    score, msa = align_guide_tree(mst, alphabet)

    # 4. Return score

    return score, msa

def sp_cost(msa, delta):
  k = len(msa)
  n = len(msa[0])

  sum = 0
  for col in range(0, n):
    characters = []
    for row in range(k):
      characters.append(msa[row][col])
    for p in range(k):
      for q in range(p+1, k):
        sum += delta[characters[p]][characters[q]]
  return sum

def align_guide_tree(G, alphabet):
    new_node_idx = G.number_of_nodes()
    score = None
    while (G.number_of_nodes() > 1):
      min_edge = min(G.edges(data=True), key=lambda x: x[2]['weight'])
      i, j, _ = min_edge
      # print(f"The minimum-weight edge is: {min_edge}")
      _, new_node_alignment = progressive_alignment(G.nodes[i]['alignment'], G.nodes[j]['alignment'], alphabet, delta_for_max)
      score = sp_cost(new_node_alignment, delta_for_min)
      # print(score, new_node_alignment)
      G.add_nodes_from([ (new_node_idx, {"alignment": new_node_alignment}) ])

      # print(f'i: {i}, j: {j}, New: {new_node_idx}')
      for neighbor in G.neighbors(i):
         if neighbor != j:
            G.add_edge(new_node_idx, neighbor, weight=G[i][neighbor]['weight'])
            # print(f'Adding edge: ({new_node_idx}, {neighbor})')

      for neighbor in G.neighbors(j):
         if neighbor != i:
            G.add_edge(new_node_idx, neighbor, weight=G[j][neighbor]['weight'])
            # print(f'Adding edge: ({new_node_idx}, {neighbor})')

      G.remove_node(i)
      G.remove_node(j)
      # print(f'Removed nodes: {i}, {j}')

      new_node_idx += 1

      # Visualize current graph
      # pos = nx.spring_layout(G)
      # nx.draw(G, pos, with_labels=True, edge_color='red')
      # labels = nx.get_edge_attributes(G, 'weight')
      # nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
      # plt.title('Current Tree')
      # plt.show()
    return score, G.nodes[new_node_idx - 1]['alignment']

from utils import get_delta_for_max, get_delta_for_min
from pairwise_alignment import pairwise_align
from itertools import product
import numpy as np
import heapq

# Global variables
data = None
k = None
sequences_length = None
heap = None
cost_table = None
delta_for_max = None
delta_for_min = None

# Return: (a, a, ..., a) to (b, b, ..., b), each of length k
def binary_iterate(k, a, b):
  return product(range(a, b + 1), repeat=k)

# Return: ( start ) to ( end )
def tuple_range(start, end):
  if (len(start) != len(end)):
    raise RuntimeError("tuple_range(): dimension of start-end don't match")

  return product(*(range(s, e + 1) for s, e in zip(start, end)))

# node: int - tuple of length k
# case: (-1 or 0, ...) - tuple of length k
def sp_cost(node, case):
  sum = 0
  characters = [data[i][node[i] - 1] if (case[i] == -1) else '-' for i in range(k)]
  for p in range(k):
    for q in range(p+1, k):
      sum += delta_for_min[characters[p]][characters[q]]
  return sum

# node: int - tuple of length k
# case: (-1 or 0, ...) - tuple of length k
def add_to_heap(old_node_tuple, case):
  new_node_np = np.subtract(old_node_tuple, case)
  new_node_tuple = tuple(new_node_np)

  # Out of index in cost_table
  if (np.any(new_node_np < 0) or np.any(new_node_np > sequences_length)):
      return

  # Calculate the cost, and add to heap
  cost = cost_table[old_node_tuple] + sp_cost(new_node_tuple, case)
  global heap
  heapq.heappush(heap, (cost, new_node_tuple, case))
  return

def carillo_lipman(raw_data, input_z):
  global data, k, sequences_length, heap, cost_table, delta_for_max, delta_for_min
  if delta_for_max == None:
      delta_for_max = get_delta_for_max()
  if delta_for_min == None:
      delta_for_min = get_delta_for_min()
  data = raw_data
  z = input_z
  k = len(data) # no. of sequences
  sequences_length = np.array([len(row) for row in data]) # length of each sequence
  cost_table = np.full(shape=(sequences_length + 1), fill_value=None) # Add 'didn't consume anything' entry at the beginning of each string
  traceback_table = np.full(shape=(cost_table.shape), fill_value=None)

  start_node = (0,) * k
  end_node = tuple(sequences_length)
  heap = [(0, start_node, None)] # Heap: (cost, node, ancestor_node)
  heapq.heapify(heap)
  while len(heap) > 0:
    # Update cost_table, traceback_table
    cost, node, min_case = heapq.heappop(heap)
    if (cost_table[node] != None) and (cost_table[node] < cost): # Already have less value for that node
      continue
    cost_table[node] = cost
    traceback_table[node] = min_case

    # Prune (if necessary)
    pairwise_sum = 0
    for i in range(k):
      for j in range(i+1, k):
        pairwise_sum += (-1) * pairwise_align(data[i][node[i]:], data[j][node[j]:], delta_for_max)[0] # Pairwise alignment of SUFFIX
    if (z < cost + pairwise_sum):
      continue # Suboptimal, do not explore from this node - pruned

    # List all cases in the recurrence
    cases = binary_iterate(k, -1, 0) # Have to create separate iterator independently, do not move outside
    all_zero_case = (0,) * k

    # Forward recurrence
    for case in cases:
      if (case == all_zero_case): # Recurrence doesn't consider this case
        continue
      add_to_heap(node, case)

  # Safeguard: Check if we reach the end_node as desired
  if (cost_table[end_node] == None):
    raise RuntimeError("Error: Cannot reach the end_node")

  # Calculate space
  total_space = np.prod(sequences_length + 1)
  pruned_space = total_space - (cost_table != None).sum()

  # Traceback and construct MSA
  current_node = end_node
  MSA_transpose = []
  while (current_node != start_node):
    min_case = traceback_table[current_node]

    characters = [data[i][current_node[i] - 1] if (min_case[i] == -1) else '-' for i in range(k)]
    MSA_transpose.insert(0, characters)

    current_node = tuple(np.add(current_node, min_case))
  MSA = np.transpose(np.array(MSA_transpose))

  return MSA, cost_table[end_node], pruned_space, total_space

def print_result(output):
  MSA, cost, pruned_space, total_space = output
  print(f'MSA:\n {MSA}')
  print(f'Cost: {cost}')
  print(f'Pruned spaces: {pruned_space}/{total_space} = {pruned_space/total_space * 100}%')

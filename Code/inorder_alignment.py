'''
Inorder Alignment:
Combine the first sequence with the second sequence, then combine the resulting MSA with the third sequence, then combine the resulting MSA with the fourth sequence, and so on
    - Another heuristic to use as z-value for Carillo-Lipman
    - Made-up algorithm (doesn't exist in real life)
'''

from progressive_alignment import align_profile_profile
from utils import get_delta_for_min

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

def inorder_alignment(alignments_input, alphabet, delta):
    """
    :param: alignments is a list of strings representing the sequences to be aligned
            Note: This is because we need to represent our single sequences as multiple alignments
            ,and multiple alignments are lists of strings
    :param: alphabet is the alphabet from which the sequences are derived
    :param: delta is a scoring function. delta(x,y) gives us the score for aligning
            character x with character y in our alphabet

    :returns: the greedy optimal multiple sequence alignment for a given set of sequences, and the score for that alignment
    """
    # Converting each sequence into MSA-like format
    alignments = []
    for seq in alignments_input:
       alignments.append([seq])

    while True:
        # Base case (When to exit the loop?) when lengtg == 1
        # YOUR CODE HERE
        if len(alignments) == 1:
          score = sp_cost(alignments[0], get_delta_for_min())
          return score, alignments[0]


        # Data structures for this iteration
        best_score = -float("inf")
        best_alignment = None
        best_m = -1
        best_n = -1

        # Combine sequence/MSA to sequence/MSA
        
        _, result_alignment = align_profile_profile(alignments[0], alignments[1], alphabet, delta)

        # Populate the list of alignments to use for the next iteration
        del alignments[0]
        del alignments[0] # Former index 1
        alignments.insert(0, result_alignment)

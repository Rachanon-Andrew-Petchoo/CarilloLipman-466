## Greedy Progressive Alignement
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

def greedy_progressive_align(alignments_input, alphabet, delta):
    """
    :param: alignments is a list of strings representing the sequences to be aligned
            Note: This is because we need to represent our single sequences as multiple alignments
            ,and multiple alignments are lists of strings
    :param: alphabet is the alphabet from which the sequences are derived
    :param: delta is a scoring function. delta(x,y) gives us the score for aligning
            character x with character y in our alphabet

    :returns: the greedy optimal multiple sequence alignment for a given set of sequences, and the score for that alignment
    """
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

        # Compute pairwise distances
        for m in range(len(alignments)):
            for n in range(m):
                # YOUR CODE HERE
                score, aligning = align_profile_profile(alignments[m], alignments[n], alphabet, delta)
                if best_score < score:
                  #maximise
                  best_score = score
                  best_alignment = aligning
                  best_n = n
                  best_m = m

        # Populate the list of alignments to use for the next iteration
        next_alignments = [best_alignment]
        for i in range(len(alignments)):
            if i!=best_m and i!=best_n:
                next_alignments.append(alignments[i])
        alignments = next_alignments


# Student test case (You may add more)

# alphabet = ['A', 'C', 'G', 'T', '-']
# delta_dna = {}
# for i in range(len(alphabet)):
#     delta_dna[alphabet[i]] = {k : v for (k,v)
#                           in zip(alphabet, [1 if alphabet[i] == alphabet[j]  else -1
#                                   for j in range(len(alphabet))]
#                          )}


# alphabet = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
# delta_protein = {}
# for i in range(len(alphabet)):
#     delta_protein[alphabet[i]] = {k : v for (k,v)
#                           in zip(alphabet, [1 if alphabet[i] == alphabet[j]  else -1
#                                   for j in range(len(alphabet))]
#                          )}


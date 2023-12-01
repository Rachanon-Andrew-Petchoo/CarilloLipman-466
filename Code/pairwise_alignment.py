UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

def traceback_global(v, w, pointers):
    '''
    Private helper function
    Used for back trace of dynamic programing table using backpointers.
    '''
    i,j = len(v), len(w)
    new_v = []
    new_w = []
    while True:
        di, dj = pointers[i][j]
        if (di,dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j-1])
        elif (di,dj) == UP:
            new_v.append(v[i-1])
            new_w.append('-')
        elif (di,dj) == TOPLEFT:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    return [''.join(new_v[::-1]),''.join(new_w[::-1])]

def max_of_three(leftValue,upValue,topleftValue):
    '''
    Private helper function
    returns: contanst TOPLEFT, LEFT, or UP
    Assumes that items is the values from a completed table
    '''

    neighbors = [leftValue, upValue, topleftValue]
    # Use max with a key function that compares the first element of each tuple
    max_location = neighbors.index(max(neighbors))
    if max_location == 0:
      return LEFT
    elif max_location == 1:
      return UP
    else:
      return TOPLEFT

def pairwise_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment,
    as a list of two strings, as computed by traceback_global. The is a global alignment of the two sequences.

    :param: v is a sequence
    :param: w is a second sequence
    :param: delta is the scoring function
    """
    M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]
    score, alignment = None, None
    # YOUR CODE HERE
 # initialize the origin with 0 score
    M[0][0] = 0;
    pointers[0][0] = ORIGIN

    # Initialize the first row and column with gap penalties
    for i in range(1, len(v)+1):
        M[i][0] = M[i-1][0] + delta[v[i-1]]['-']
        pointers[i][0] = UP
    for j in range(1, len(w)+1):
        M[0][j] = M[0][j-1] + delta['-'][w[j-1]]
        pointers[0][j] = LEFT

    # fill in remaining table starting a col and row 1. (note: for v and w, subtract 1 from i, j b/c they don't have leading zero)
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
          score_TL   = M[i-1][j-1] + delta[v[i-1]][w[j-1]]   # match or mismatch given that delta will score
          score_UP   = M[i-1][j]   + delta[v[i-1]]['-']      # deletion
          score_LEFT = M[i][j-1]   + delta['-'][w[j-1]]      # insertion
          max_score = max(score_LEFT, score_UP, score_TL)
          M[i][j] = max_score
          pointers[i][j] = max_of_three(score_LEFT, score_UP, score_TL)

    score = M[len(v)][len(w)]

    alignment = traceback_global(v,w, pointers)
    return score, alignment

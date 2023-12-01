UP = (-1, 0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)

def compute_profile(alignment, alphabet):
    """
    Given an alphabet an a multiple sequence alignment in that alphabet,
    computes and returns its profile representation

    :param: alignment is a list of lists of characters in the alphabet
    :param: alphabet is a list of characters in the alphabet from which the strings are
            constructed
    :return: a dictionary where dict[x][i] is the frequency of the character
             x in the i-th position of the alignment.
    """

    if not alignment:
        return {}

    n = len(alignment)
    length = len(alignment[0])
    profile = {}

    # init profile
    for letter in alphabet:
      profile[letter] = []

    # for each column in MSA
    for i in range(length):
      letters = [sequence[i] for sequence in alignment] # create a list of letters in i column
      letters = "".join(letters)                        # convert list to a string
      for letter in alphabet:
        count = letters.count(letter)                   # get number of letter
        profile[letter].append(count / n)               # store frequency of letter in profile
       # print(f"freq={profile[letter]} for {letter}")

    return profile

def compute_tau(profile, alphabet, delta):
    """
    Given a profile, an alphabet and a scoring function for that alphabet,
    returns the scoring function for aligning a character in the alphabet
    to a column in the profile

    :param: profile is the profile representation of the multiple sequence
            we are aligning against
    :param: alphabet is the alphabet of characters that compose our sequences
    :param: delta is the scoring function between characters in our alphabet

    :return: The scoring function tau such that tau[x][i] is the score for aligning
             character x with column i of the profile.
    """
    tau = {}

    for letter in alphabet:
      tau[letter] = []

    length = len(profile[alphabet[0]])     # num of columns in the MSA
    # size_of_alphabet = len(profile)   # probably could just do len(alphabet) but did this way incase profile not same

    for j in range(length):
      for letter in alphabet:
        sum = 0
        for key in profile.keys():    # i is the index of the column of the MSA
          # print(f"i={i}, key={key} pro={profile}")
          sum += profile[key][j] * delta[letter][key]
        tau[letter].append(sum)
    return tau

def compute_sigma(p,q,alphabet, delta):
    """
    :param: p is the profile for the first multiple alignment
    :param: q is the profile for the second multiple alignment
    :param: alphabet is the list of all characters in our sequences
    :param: delta is the scoring function for aligning characters in our alphabet

    :returns: a list of lists sigma such that sigma[i][j] is the score for aligning column
              i of p with column j of q
    """
    sigma = []
    size_p = len(p[alphabet[0]])
    size_q = len(q[alphabet[0]])

    # Loop through the rows and columns of the sigma matrix
    for i in range(size_p):
        row = []
        for j in range(size_q):
            element = 0
            # Loop through the alphabet and calculate the element
            for x in alphabet:
                for y in alphabet:
                    element += p[x][i] * q[y][j] * delta[x][y]
            # Append the element to the row
            row.append(element)
        # Append the row to the sigma matrix
        sigma.append(row)

    return sigma

def traceback(aln1, aln2, pointers):
    i = len(aln1[0])-1
    j = len(aln2[0])-1
    new_al1 = [list(v) for v in aln1]
    new_al2 = [list(w) for w in aln2]
    while True:
        di, dj = pointers[i][j]
        if (di, dj) == LEFT:
            for seq1 in new_al1:
                seq1.insert(i, '-')
        if (di, dj) == UP:
            for seq2 in new_al2:
                seq2.insert(j, '-')
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    new_alignment = []
    for seq in new_al1:
        new_alignment.append(''.join(seq))
    for seq in new_al2:
        new_alignment.append(''.join(seq))
    return new_alignment

def align_sequence_profile(alignment, sequence, alphabet, delta):
    """
    This function aligns a sequence against a multiple sequence alignment

    :param: alignment is the multiple sequence alignment are aligning against.
            This is a list of list of characters
    :param: sequence is the new sequence we are aligning to the multiple alignment.
            This is a list of characters
    :param: alphabet is a list of characters that could compose the sequences in
            the alignments.
    :param: delta is the scoring function for aligning characters in our alphabet.
            delta[x][y] is the score for aligning the characters x and y.


    :return: a list of lists of characters in the alphabet, representing the
             new multiple sequence alignment
    """
    # Base case when there is an empty multiple alignment
    if not alignment:
        return [sequence]
    M = [[0 for _ in range(len(alignment[0]))] for _ in range(len(sequence))]
    pointers = [[(0,0) for _ in range(len(alignment[0]))] for _ in range(len(sequence))]
    score = None

    profile = compute_profile(alignment, alphabet)
    tau = compute_tau(profile, alphabet, delta)

    for i in range(len(sequence)):
        for j in range(len(alignment[0])):
            if i == 0 and j == 0:
                M[i][j] = 0
            elif i == 0:
                M[i][j] = M[i][j-1] + tau['-'][j-1]
                pointers[i][j] = LEFT
            elif j == 0:
                sequence[i-1]
                M[i][j] = M[i-1][j] + delta[sequence[i-1]]['-']
                pointers[i][j] = UP
            else:
                best_sub = max([(LEFT, M[i][j-1] + tau['-'][j-1]),
                               (UP, M[i-1][j] + delta[sequence[i-1]]['-']),
                               (TOPLEFT, M[i-1][j-1] + tau[sequence[i-1]][j-1])], key = lambda x: x[1])
                pointers[i][j] = best_sub[0]
                M[i][j] = best_sub[1]

    score = M[-1][-1]
    return score, traceback([sequence], alignment, pointers)

def align_profile_profile(aln1, aln2, alphabet, delta):
    """
    :param: aln1 is a list of lists representing the first multiple alignment
    :param: aln2 is a list of lists representing the second multiple alignment
    :param: alphabet is the alphabet from which the sequences are derived
    :param: delta is a scoring function. delta(x,y) gives us the score for aligning
            character x with character y in our alphabet

    :returns: the optimal score and the optimal multiple alignment for the two input alignments.
    """
    # Base case when there is an empty multiple alignment
    if not aln1 and not aln2:
        return []
    elif not aln2:
        return aln1
    elif not aln1:
        return aln2

    S = [[0 for j in range(len(aln2[0])+1)] for i in range(len(aln1[0])+1)]
    pointers = [[(0,0) for j in range(len(aln2[0])+1)] for i in range(len(aln1[0])+1)]
    score = None

    # Compute profiles and scoring functions here
    p = compute_profile(aln1, alphabet)
    q = compute_profile(aln2, alphabet)
    tau1 = compute_tau(p, alphabet, delta)
    tau2 = compute_tau(q, alphabet, delta)
    sigma = compute_sigma(p,q,alphabet, delta)

    for i in range(1, len(aln1[0]) + 1):
      S[i][0] = S[i - 1][0] + tau1['-'][i - 1]
      pointers[i][0] = UP
    for i in range (1, len(aln2[0]) + 1):
      S[0][i] = S[0][i - 1] + tau2['-'][i - 1]
      pointers[0][i] = LEFT

    ## table from slides
    for i in range(1, len(aln1[0]) + 1):
      for j in range(1, len(aln2[0]) + 1):
        l = S[i][j - 1] + tau2['-'][j - 1]
        u = S[i - 1][j] + tau1['-'][i - 1]
        d = S[i - 1][j - 1] + sigma[i - 1][j - 1]
        best_score = max(l, u, d)
        if (best_score == l):
          pointers[i][j] = LEFT
        elif (best_score == u):
          pointers[i][j] = UP
        else:
          pointers[i][j] = TOPLEFT

        S[i][j] = best_score

    score = S[-1][-1]
    return score, traceback(aln1, aln2, pointers)

def progressive_alignment(aln1, aln2, alphabet, delta):
    """
    :param: aln1 is a list of lists representing the first multiple alignment
    :param: aln2 is a list of lists representing the second multiple alignment
    :param: alphabet is the alphabet from which the sequences are derived
    :param: delta is a scoring function. delta(x,y) gives us the score for aligning
            character x with character y in our alphabet

    :returns: the optimal multiple alignment for the two input alignments.
    """
    # Base case when there is an empty multiple alignment
    if not aln1 and not aln2:
        return []
    elif not aln2:
        return aln1
    elif not aln1:
        return aln2

    S = [[0 for j in range(len(aln2[0])+1)] for i in range(len(aln1[0])+1)]
    pointers = [[(0,0) for j in range(len(aln2[0])+1)] for i in range(len(aln1[0])+1)]
    score = None

    # Compute profiles and scoring functions here
    p = compute_profile(aln1, alphabet)
    q = compute_profile(aln2, alphabet)
    tau1 = compute_tau(p, alphabet, delta)
    tau2 = compute_tau(q, alphabet, delta)
    sigma = compute_sigma(p, q, alphabet, delta)

    for i in range(len(aln1[0])+1):
        for j in range(len(aln2[0])+1):
            if i == 0 and j == 0:
                S[i][j] = 0
            elif i == 0:
                S[i][j] = S[i][j-1] + tau2['-'][j-1]
                pointers[i][j] = LEFT
            elif j == 0:
                S[i][j] = S[i-1][j] + tau1['-'][i-1]
                pointers[i][j] = UP
            else:
                best_sub = max([(LEFT,S[i][j-1] + tau2['-'][j-1]) ,
                               (UP, S[i-1][j] + tau1['-'][i-1]),
                               (TOPLEFT, S[i-1][j-1] + sigma[i-1][j-1])], key = lambda x: x[1])
                pointers[i][j] = best_sub[0]
                S[i][j] = best_sub[1]

    score = S[-1][-1]
    return score, traceback(aln1, aln2, pointers)

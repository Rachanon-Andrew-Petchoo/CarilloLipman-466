import random

# Global variables
keys = None
delta_for_min = None
delta_for_max = None

def get_DNA_keys():
   global keys
   if keys != None:
      return keys
   
   keys = ['A', 'C', 'T', 'G', '-']
   return keys

# Score function for functions that want to minimize
def get_delta_for_min():
    global keys, delta_for_min
    if delta_for_min != None:
        return delta_for_min
    
    if keys == None:
        get_DNA_keys()
    
    delta_for_min = {}
    for i in range(len(keys)):
        if keys[i] != '-':
            delta_for_min[keys[i]] = {k : v for (k,v) in zip(keys, [-3 if keys[i] == keys[j]  else -1 for j in range(len(keys))])}
    delta_for_min['-'] = {}
    for k in keys:
        delta_for_min[k]['-'] = 1
        delta_for_min['-'][k] = 1
    delta_for_min['-']['-'] = 0

    return delta_for_min

# Score function for functions that want to maximize (opposite sign with get_delta_for_min())

def get_delta_for_max():
    global delta_for_min, delta_for_max
    if delta_for_max != None:
        return delta_for_max
    
    if delta_for_min == None:
        get_delta_for_min()
    
    delta_for_max = {}
    for k in delta_for_min:
      delta_for_max[k] = {}
      for v in delta_for_min[k]:
        delta_for_max[k][v] = - delta_for_min[k][v]
    
    return delta_for_max

# Function to generate random DNA sequences
def generate_random_test_dna(seq_len, num_seq):
    nucleotides = ['A', 'C', 'G', 'T']
    return [''.join(random.choices(nucleotides, k=seq_len)) for i in range(num_seq)]

# Function to generate random Protein sequences
def generate_random_test_protein(seq_len, num_seq):
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    return [''.join(random.choices(amino_acids, k=seq_len)) for i in range(num_seq)]
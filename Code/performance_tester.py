from carillo_lipman import carillo_lipman
from greedy_progressive_alignment import greedy_progressive_align
from feng_and_doolittle import feng_doolittle_alignment
from inorder_alignment import inorder_alignment
from utils import generate_random_test_dna, get_delta_for_max
import time
import matplotlib.pyplot as plt
import numpy as np

# Track runtime for DNA data
runtime_dna_feng_doolittle = []
runtime_dna_greedy = []
runtime_dna_inorder = []

# Track space for DNA data
space_dna_feng_doolittle = []
space_dna_greedy = []
space_dna_inorder = []

# Track [ [sequence length, space], [_, _], ..., [_,_] ] for DNA data
length_space_dna_feng_doolittle = []
length_space_dna_greedy = []
length_space_dna_inorder = []

# Track [ [sequence length, runtime], [_, _], ..., [_,_] ] for DNA data
length_runtime_dna_feng_doolittle = []
length_runtime_dna_greedy = []
length_runtime_dna_inorder = []

def dna_tester(iters, num_seq, seq_len, alphabet, delta):
    print(f'\nDNA Tester: {num_seq} sequences with length {seq_len}')
    global runtime_dna_feng_doolittle, runtime_dna_greedy, runtime_dna_inorder
    global space_dna_feng_doolittle, space_dna_greedy, space_dna_inorder
    global length_space_dna_feng_doolittle, length_space_dna_greedy, length_space_dna_inorder
    global length_runtime_dna_feng_doolittle, length_runtime_dna_greedy, length_runtime_dna_inorder

    for i in range(iters):
        test_case = generate_random_test_dna(seq_len, num_seq)

        # Greedy:
        score, _ = greedy_progressive_align(test_case, alphabet, delta)
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()
        runtime = end - start
        space_percentage = (1 - (pruned_space/total_space)) * 100

        runtime_dna_greedy.append(runtime)
        space_dna_greedy.append(space_percentage)
        length_space_dna_greedy.append([seq_len, space_percentage])
        length_runtime_dna_greedy.append([seq_len, runtime])

        print(f'Greedy score: {score}, {space_percentage}%')

        # Feng and Doolittle:
        score, _ = feng_doolittle_alignment(test_case, alphabet)
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()
        runtime = end - start
        space_percentage = (1 - (pruned_space/total_space)) * 100

        runtime_dna_feng_doolittle.append(runtime)
        space_dna_feng_doolittle.append(space_percentage)
        length_space_dna_feng_doolittle.append([seq_len, space_percentage])
        length_runtime_dna_feng_doolittle.append([seq_len, runtime])

        print(f'Feng and Doolittle score: {score}, {space_percentage}%')

        # Inorder:
        score, _ = inorder_alignment(test_case, alphabet, delta)
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()
        runtime = end - start
        space_percentage = (1 - (pruned_space/total_space)) * 100

        runtime_dna_inorder.append(runtime)
        space_dna_inorder.append(space_percentage)
        length_space_dna_inorder.append([seq_len, space_percentage])
        length_runtime_dna_inorder.append([seq_len, runtime])

        print(f'Inorder score: {score}, {space_percentage}%')

def performance_test():
    global runtime_dna_feng_doolittle, runtime_dna_greedy, runtime_dna_inorder
    global space_dna_feng_doolittle, space_dna_greedy, space_dna_inorder

    alphabet_dna = ['A', 'C', 'G', 'T', '-']
    delta_dna = get_delta_for_max()

    # Create a plot for Space vs. Runtime in DNA
    # NOTE: DO NOT SWAP ORDER OF CODE - only 1 setting of 'sequence length' and 'no. of sequences' should be populated
    dna_tester(20, 5, 9, alphabet_dna, delta_dna)
    concat_runtime_dna = runtime_dna_feng_doolittle + runtime_dna_greedy + runtime_dna_inorder
    concat_space_dna = space_dna_feng_doolittle + space_dna_greedy + space_dna_inorder
    plt.scatter(concat_space_dna, concat_runtime_dna)
    plt.title('Space vs. Runtime')
    plt.xlabel('Space used (percentage)')
    plt.ylabel('Runtime')
    plt.show()
    
    # iters, no. of sequences, sequence length
    # NOTE: DO NOT SWAP ORDER OF CODE - 'sequence length' vs. 'runtime / space' plots depends on it
    dna_tester(10, 5, 5, alphabet_dna, delta_dna)
    dna_tester(10, 5, 6, alphabet_dna, delta_dna)
    dna_tester(10, 5, 7, alphabet_dna, delta_dna)
    dna_tester(20, 5, 8, alphabet_dna, delta_dna)

    # Create a plot for DNA sequence length vs. space
    fig, axs = plt.subplots(2, 2)
    # Greedy
    length_space_dna_greedy_np = np.array(length_space_dna_greedy, dtype=int)
    axs[0, 0].scatter(length_space_dna_greedy_np[:,0], length_space_dna_greedy_np[:,1])
    axs[0, 0].set_title('Greedy Progressive')
    # Feng and Doolittle
    length_space_dna_feng_doolittle_np = np.array(length_space_dna_feng_doolittle)
    axs[0, 1].scatter(length_space_dna_feng_doolittle_np[:,0], length_space_dna_feng_doolittle_np[:,1])
    axs[0, 1].set_title('Feng & Doolittle')
    # Inorder
    length_space_dna_inorder_np = np.array(length_space_dna_inorder)
    axs[1, 0].scatter(length_space_dna_inorder_np[:,0], length_space_dna_inorder_np[:,1])
    axs[1, 0].set_title('Inorder alignment')
    # Combine plots
    for ax in axs.flat:
        ax.set(xlabel='Sequence length', ylabel='Space used (percentage)')
    for ax in axs.flat:
        ax.label_outer()
    plt.show()

    # Create a plot for DNA sequence length vs. runtime
    fig, axs = plt.subplots(2, 2)
    # Greedy
    length_runtime_dna_greedy_np = np.array(length_runtime_dna_greedy, dtype=int)
    axs[0, 0].scatter(length_runtime_dna_greedy_np[:,0], length_runtime_dna_greedy_np[:,1])
    axs[0, 0].set_title('Greedy Progressive')
    # Feng and Doolittle
    length_runtime_dna_feng_doolittle_np = np.array(length_runtime_dna_feng_doolittle)
    axs[0, 1].scatter(length_runtime_dna_feng_doolittle_np[:,0], length_runtime_dna_feng_doolittle_np[:,1])
    axs[0, 1].set_title('Feng & Doolittle')
    # Inorder
    length_runtime_dna_inorder_np = np.array(length_runtime_dna_inorder)
    axs[1, 0].scatter(length_runtime_dna_inorder_np[:,0], length_runtime_dna_inorder_np[:,1])
    axs[1, 0].set_title('Inorder alignment')
    # Combine plots
    for ax in axs.flat:
        ax.set(xlabel='Sequence length', ylabel='Runtime')
    for ax in axs.flat:
        ax.label_outer()
    plt.show()

    # NOTE: DO NOT SWAP ORDER OF CODE - 'sequence length' vs. 'runtime / space' plots depends on it
    dna_tester(10, 6, 5, alphabet_dna, delta_dna)
    dna_tester(20, 6, 6, alphabet_dna, delta_dna)
    dna_tester(20, 6, 7, alphabet_dna, delta_dna)
    dna_tester(20, 6, 8, alphabet_dna, delta_dna)
    dna_tester(20, 7, 5, alphabet_dna, delta_dna)

    # Create a box plot for DNA space
    space_dna = [space_dna_feng_doolittle,
                 space_dna_greedy,
                 space_dna_inorder]
    
    plt.boxplot(space_dna)
    plt.title('Space used for each heuristic (DNA data)')
    plt.xlabel('Heuristics')
    plt.ylabel('Space used (percentage)')
    plt.xticks([1, 2, 3], ['Feng & Doolittle', 'Greedy Progressive', 'Inorder'])
    plt.show()

    # Create a box plot for DNA runtime
    runtime_dna = [runtime_dna_feng_doolittle, 
                   runtime_dna_greedy,
                   runtime_dna_inorder]

    plt.boxplot(runtime_dna)
    plt.title('Runtime for each heuristic (DNA data)')
    plt.xlabel('Heuristics')
    plt.ylabel('Runtime')
    plt.xticks([1, 2, 3], ['Feng & Doolittle', 'Greedy Progressive', 'Inorder'])
    plt.show()

if __name__ == "__main__":
    performance_test()
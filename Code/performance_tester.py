from carillo_lipman import carillo_lipman
from greedy_progressive_alignment import greedy_progressive_align
from feng_and_doolittle import feng_doolittle_alignment
from utils import generate_random_test_dna, generate_random_test_protein, get_delta_for_max
import time
import matplotlib.pyplot as plt

# Track runtime for DNA data
runtime_dna_feng_doolittle = []
runtime_dna_greedy = []
runtime_dna_ClustalW = []

# Track space for DNA data
space_dna_feng_doolittle = []
space_dna_greedy = []
space_dna_ClustalW = []

# Track runtime for Protein data
runtime_protein_feng_doolittle = []
runtime_protein_greedy = []
runtime_protein_ClustalW = []

# Track space for Protein data
space_protein_feng_doolittle = []
space_protein_greedy = []
space_protein_ClustalW = []

def dna_tester(iters, seq_len, num_seq, alphabet, delta):
    global runtime_dna_feng_doolittle, runtime_dna_greedy, runtime_dna_ClustalW
    global space_dna_feng_doolittle, space_dna_greedy, space_dna_ClustalW

    for i in range(iters):
        test_case = generate_random_test_dna(seq_len, num_seq)
        # Greedy:
        score, _ = greedy_progressive_align(test_case, alphabet, delta)
        print(f'Greedy score: {score}')
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()

        runtime_dna_greedy.append(end - start)
        space_dna_greedy.append((1 - (pruned_space/total_space))*100)

        # Feng and Doolittle:
        score, _ = feng_doolittle_alignment(test_case, alphabet)
        print(f'Feng and Doolittle score: {score}')
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()

        runtime_dna_feng_doolittle.append(end - start)
        space_dna_feng_doolittle.append((1 - (pruned_space/total_space))*100)

        # ClustalW:
        # TODO: NICK PUT IT HERE in the same format

def protein_tester(iters, seq_len, num_seq, alphabet, delta):
    global runtime_protein_feng_doolittle, runtime_protein_greedy, runtime_protein_ClustalW
    global space_protein_feng_doolittle, space_protein_greedy, space_protein_ClustalW

    for i in range(iters):
        test_case = generate_random_test_protein(seq_len, num_seq)

        # Greedy:
        score, _ = greedy_progressive_align(test_case, alphabet, delta)
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()

        runtime_protein_greedy.append(end - start)
        space_protein_greedy.append((1 - (pruned_space/total_space))*100)

        # Feng and Doolittle:
        score, _ = feng_doolittle_alignment(test_case, alphabet)
        start = time.time()
        _, _, pruned_space, total_space = carillo_lipman(test_case, score)
        end = time.time()

        runtime_protein_feng_doolittle.append(end - start)
        space_protein_feng_doolittle.append((1 - (pruned_space/total_space))*100)

        # ClustalW:
        # TODO: NICK PUT IT HERE in the same format

def performance_test():
    global runtime_dna_feng_doolittle, runtime_dna_greedy, runtime_dna_ClustalW
    global space_dna_feng_doolittle, space_dna_greedy, space_dna_ClustalW
    global runtime_protein_feng_doolittle, runtime_protein_greedy, runtime_protein_ClustalW
    global space_protein_feng_doolittle, space_protein_greedy, space_protein_ClustalW

    alphabet_dna = ['A', 'C', 'G', 'T', '-']
    delta_dna = get_delta_for_max()
    
    # iters, sequence length, no. of sequences
    dna_tester(10, 5, 5, alphabet_dna, delta_dna)
    dna_tester(10, 6, 6, alphabet_dna, delta_dna)
    # dna_tester(10, 10, 5, alphabet_dna, delta_dna)
    # dna_tester(10, 15, 5, alphabet_dna, delta_dna)
    # dna_tester(10, 5, 10, alphabet_dna, delta_dna)
    # dna_tester(10, 10, 10, alphabet_dna, delta_dna)
    # dna_tester(10, 20, 10, alphabet_dna, delta_dna)

    #TODO: Use other delta for Protein
    alphabet_protein = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
    delta_protein = {}
    for i in range(len(alphabet_protein)):
        delta_protein[alphabet_protein[i]] = {k : v for (k,v)
                            in zip(alphabet_protein, [1 if alphabet_protein[i] == alphabet_protein[j]  else -1
                                    for j in range(len(alphabet_protein))]
                            )}

    # protein_tester(10, 20, 5, alphabet_protein, delta_protein)
    # protein_tester(10, 20, 10, alphabet_protein, delta_protein)
    # protein_tester(10, 20, 15, alphabet_protein, delta_protein)
    # protein_tester(10, 30, 5, alphabet_protein, delta_protein)
    # protein_tester(10, 30, 10, alphabet_protein, delta_protein)
    # protein_tester(10, 30, 15, alphabet_protein, delta_protein)

    # Create a box plot for DNA runtime
    runtime_dna = [runtime_dna_feng_doolittle, 
                   runtime_dna_greedy, 
                   runtime_dna_ClustalW]

    plt.boxplot(runtime_dna)
    plt.title('Runtime for each heuristic (DNA data)')
    plt.xlabel('Heuristics')
    plt.ylabel('Runtime')
    plt.xticks([1, 2, 3], ['Feng and Doolittle', 'Greedy Progressive', 'ClustalW'])
    plt.show()

    # Create a box plot for DNA space
    space_dna = [space_dna_feng_doolittle, 
                   space_dna_greedy, 
                   space_dna_ClustalW]
    
    plt.boxplot(space_dna)
    plt.title('Space used for each heuristic (DNA data)')
    plt.xlabel('Heuristics')
    plt.ylabel('Space used (percentage)')
    plt.xticks([1, 2, 3], ['Feng and Doolittle', 'Greedy Progressive', 'ClustalW'])
    plt.show()

    # # Create a box plot for Protein runtime
    # runtime_protein = [runtime_protein_feng_doolittle,
    #                    runtime_protein_greedy,
    #                    runtime_protein_ClustalW]

    # plt.boxplot(runtime_protein)
    # plt.title('Runtime for each heuristic (Protein data)')
    # plt.xlabel('Heuristics')
    # plt.ylabel('Runtime')
    # plt.xticks([1, 2, 3], ['Feng and Doolittle', 'Greedy Progressive', 'ClustalW'])
    # plt.show()

    # # Create a box plot for DNA space
    # space_protein = [space_protein_feng_doolittle,
    #                  space_protein_greedy,
    #                  space_protein_ClustalW]
    
    # plt.boxplot(space_protein)
    # plt.title('Space used for each heuristic (Protein data)')
    # plt.xlabel('Heuristics')
    # plt.ylabel('Space used (percentage)')
    # plt.xticks([1, 2, 3], ['Feng and Doolittle', 'Greedy Progressive', 'ClustalW'])
    # plt.show()

    ############################################################################
    # # take averages?
    # results_final2 = [space_dna_feng_doolittle, space_dna_greedy, space_dna_ClustalW, prune_protein_fd, prune_protein_greedy, prune_protein_ClustalW]


    # # Create a scatter plot
    # plt.boxplot(results_final2)

    # # Customize the plot if needed
    # plt.title('Box Plot of DNA and Protein Runtimes')
    # plt.xlabel('Heuristics')
    # plt.ylabel('Runtimes')
    # plt.xticks([1, 2, 3], ['Feng and Doolittle', 'Greedy Progressive', 'ClustalW'])

    # # Show the plot
    # plt.show()

    # TODO: still:
    # get box plots to work
    # then worry about taking averages and generating scatter plots.
    # get scatter plot to work for pruned space vs runtime
    ############################################################################

if __name__ == "__main__":
    performance_test()
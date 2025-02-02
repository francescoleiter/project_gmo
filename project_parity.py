
import pyRA

import pygad
import numpy as np

best_fitness = []

N_REACTIONS = 5

num_generations = 500
sol_per_pop = 200
num_parents_mating = 70

#           initial state +        5 reactions          +   2 final states
num_genes = pyRA.ALPHABET + 3*pyRA.ALPHABET*N_REACTIONS + 2*pyRA.ALPHABET

#initial values from 0 to 1
init_range_low = 0
init_range_high = 2

# GA parameters
parent_selection_type = "sss" # steady state selection
keep_parents = -1             # keep parents
keep_elitism = 0
crossover_type = "single_point" # one-point crossover

p_geom_distribution = 0.5 # probability of the geometric distribution used for generating reactions

# improperly used parameters that may vary during computation
mutation_max_val = 0.0 # unused
mutation_probability = 1.0  #used to set the number of mutation to produce
mutation_min_val = 0.1 # probability of being non-zero for reactants and products

inh_prob = 0.02 # probability of an inhibitor to be set on (up to contradictions)
probs = np.array([0.1,0.8,0.1]) # probabilities of mutating the initial state, the reactions, the final states

# generating binary sequences
def generate_sequence(k):
    if k == 0:
        return [[]]
    return [seq + [0] for seq in generate_sequence(k-1)] + [ seq + [1] for seq in generate_sequence(k-1)]

# number of bits
bits_parity = 9
# k for k-even parity problem
k_parity = 3

# initializing k-even and k-odd sequences
acc_seqs = [seq for seq in generate_sequence(bits_parity) if not (sum(seq)%k_parity)]
rej_seqs = [seq for seq in generate_sequence(bits_parity) if (sum(seq)%k_parity)]


# function to display the process while running and save fitness values
def on_gen(ga_instance):
    best_fitness.append(ga_instance.best_solution()[1])
    print("Generation : ", ga_instance.generations_completed)
    print("Fitness of the best solution :", ga_instance.best_solution()[1])
    if not (ga_instance.generations_completed % 10):
        #pyRA.FGenome(pyRA.CVect(ga_instance.best_solution()[0])).print()
        pyRA.AGenome(pyRA.CVect(ga_instance.population[ga_instance.best_solution()[2]]), N_REACTIONS).print()


# defining fitness using the library pyRA
def fitness_function_seq(ga_instance, individual, idx, acc_seqs, rej_seqs):
    samples = list(zip(acc_seqs+rej_seqs, [True]*len(acc_seqs) + [False]*len(rej_seqs)))
    words, accs = zip(*samples)
    return pyRA.eval_fitness_seq(pyRA.AGenome(pyRA.CVect(individual), N_REACTIONS),list(map(pyRA.CVect,words)), accs, 40)

# for pure RA
def fitness_function_seq_pure(ga_instance, individual, idx, acc_seqs, rej_seqs):
    samples = list(zip(acc_seqs+rej_seqs, [True]*len(acc_seqs) + [False]*len(rej_seqs)))
    words, accs = zip(*samples)
    return pyRA.eval_fitness_seq_pure(pyRA.AGenome(pyRA.CVect(individual), N_REACTIONS),list(map(pyRA.CVect,words)), accs, 40)

# mutation operator
def mutation_func_seq(offspring, ga_instance):
    for _ in range(int(1/ga_instance.mutation_probability)):
        pyRA.agenome_2Dmutation_ow(offspring, N_REACTIONS,probs, 1.0, p_geom_distribution,
                                ga_instance.random_mutation_min_val, inh_prob)
    return offspring

for index in range(10):
    # reset stored values
    best_fitness = []
    # initialize the class
    ga_problem = pygad.GA(
                           gene_type=np.uint16,
                           keep_elitism=keep_elitism,
                           num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=lambda ga_instance, individual, idx: fitness_function_seq(ga_instance, individual, idx, acc_seqs, rej_seqs),
                           sol_per_pop=sol_per_pop,
                           num_genes=num_genes,
                           init_range_low=init_range_low,
                           init_range_high=init_range_high,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           random_mutation_min_val=mutation_min_val,
                           random_mutation_max_val=mutation_max_val,
                           crossover_type=crossover_type,
                           mutation_type=mutation_func_seq,
                           mutation_probability=mutation_probability,
                           on_generation=on_gen,
                           stop_criteria=f"reach_{len(acc_seqs)+len(rej_seqs)}"
                            )

    # run the genetic algorithm
    ga_problem.run()
    # save the values stored
    best_fitness = np.array(best_fitness)
    np.save(f"{k_parity}-parity_{index}", best_fitness)

    # print results
    print("Best individual:")
    x = pyRA.AGenome(pyRA.CVect(ga_problem.best_solution()[0]), N_REACTIONS)
    x.print()
    automaton = pyRA.read_agenome(x)

    # test on accepting and rejecting sequences
    print("Testing on input:")
    for seq in acc_seqs:
        automaton = pyRA.read_agenome(x)
        print("reading:", seq)
        for s in seq:
            print("input: ", s)
            automaton.read_single_input(s)
            automaton.verbose_react()
        while not automaton.is_halted():
            automaton.verbose_react()
        print(f"---------{automaton.accepted()}---------")
    for seq in rej_seqs:
        automaton = pyRA.read_agenome(x)
        print("reading:", seq)
        for s in seq:
            print("input:", s)
            automaton.read_single_input(s)
            automaton.verbose_react()
            print()
        while not automaton.is_halted():
            automaton.verbose_react()
        print(f"---------{automaton.accepted()}---------")

    # reset stored values
    best_fitness = []
    # for pure RA
    ga_problem = pygad.GA(
                           gene_type=np.uint16,
                           keep_elitism=keep_elitism,
                           num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=lambda ga_instance, individual, idx: fitness_function_seq_pure(ga_instance, individual, idx, acc_seqs, rej_seqs),
                           sol_per_pop=sol_per_pop,
                           num_genes=num_genes,
                           init_range_low=init_range_low,
                           init_range_high=init_range_high,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           random_mutation_min_val=mutation_min_val,
                           random_mutation_max_val=mutation_max_val,
                           crossover_type=crossover_type,
                           mutation_type=mutation_func_seq,
                           mutation_probability=mutation_probability,
                           on_generation=on_gen,
                           stop_criteria=f"reach_{len(acc_seqs)+len(rej_seqs)}"
                            )

    # run the genetic algorithm
    ga_problem.run()
    # save stored values
    best_fitness = np.array(best_fitness)
    np.save(f"{k_parity}-parity_pure_{index}", best_fitness)

    # print results
    print("Best individual:")
    x = pyRA.AGenome(pyRA.CVect(ga_problem.best_solution()[0]), N_REACTIONS)
    x.print()
    automaton = pyRA.read_apgenome(x)

    # test on accepting and rejecting sequences
    print("Testing on input:")
    for seq in acc_seqs:
        automaton = pyRA.read_apgenome(x)
        print("reading:", seq)
        for s in seq:
            print("input: ", s)
            automaton.read_single_input(s)
            automaton.verbose_react()
        while not automaton.is_halted():
            automaton.verbose_react()
        print(f"---------{automaton.accepted()}---------")
    for seq in rej_seqs:
        automaton = pyRA.read_apgenome(x)
        print("reading:", seq)
        for s in seq:
            print("input:", s)
            automaton.read_single_input(s)
            automaton.verbose_react()
            print()
        while not automaton.is_halted():
            automaton.verbose_react()
        print(f"---------{automaton.accepted()}---------")


import pyRA

import pygad
import numpy as np


N_REACTIONS = 5

num_generations = 1000
sol_per_pop = 100
num_parents_mating = 70

best_fitness = np.zeros(num_generations)

#           initial state +        5 reactions         +  3 accepting states
num_genes = pyRA.ALPHABET + 3*pyRA.ALPHABET*N_REACTIONS + 3*pyRA.ALPHABET

#initial values from 0 to 1
init_range_low = 0
init_range_high = 2

# GA parameters
parent_selection_type = "rank"
keep_parents = 0
keep_elitism = 10
crossover_type = "single_point"

p_geom_distribution = 0.75

# improperly used parameters
mutation_max_val = 0.0
mutation_probability = 0.3  #used to set the number of mutation to produce
mutation_min_val = 0.12 # probability of being non-zero for reactants and products

inh_prob = 0.1 # probability of an inhibitor to be set on (up to contradictions)
probs = np.array([0.1,0.8,0.1]) # probabilities of mutating the initial state, the reactions, the final states

# generating binary sequences
def generate_sequence(k):
    if k == 0:
        return [[]]
    return [seq + [0] for seq in generate_sequence(k-1)] + [ seq + [1] for seq in generate_sequence(k-1)]

# number of bits
k_parity = 9

# initializing even and odd sequences
acc_seqs = [seq for seq in generate_sequence(k_parity) if not (sum(seq)%2)]
rej_seqs = [seq for seq in generate_sequence(k_parity) if (sum(seq)%2)]

# defining fitness using the library pyRA
def fitness_function_seq(ga_instance, individual, idx, acc_seqs, rej_seqs):
    samples = list(zip(acc_seqs+rej_seqs, [True]*len(acc_seqs) + [False]*len(rej_seqs)))
    words, accs = zip(*samples)
    return pyRA.eval_fitness_seq(pyRA.AGenome(pyRA.Vect(individual), N_REACTIONS),list(map(pyRA.Vect,words)), accs, 40)

# for pure RA
def fitness_function_seq_pure(ga_instance, individual, idx, acc_seqs, rej_seqs):
    samples = list(zip(acc_seqs+rej_seqs, [True]*len(acc_seqs) + [False]*len(rej_seqs)))
    words, accs = zip(*samples)
    return pyRA.eval_fitness_seq_pure(pyRA.AGenome(pyRA.Vect(individual), N_REACTIONS),list(map(pyRA.Vect,words)), accs, 40)

# mutation operator
def mutation_func_seq(offspring, ga_instance):
    for _ in range(int(1/ga_instance.mutation_probability)):
        pyRA.agenome_mutation_2Darr(offspring, N_REACTIONS,probs, 1.0, p_geom_distribution,
                                ga_instance.random_mutation_min_val, inh_prob)
    return offspring

# function to display the process and update the number of mutations occurring
def on_gen(ga_instance):
    if ga_instance.best_solution()[1] <= np.max(ga_instance.previous_generation_fitness):
        ga_instance.mutation_probability = ga_instance.mutation_probability*0.9 + 0.2*0.1
    else:
        ga_instance.mutation_probability = 1.0
    if not (ga_instance.generations_completed%10):
        print("Generation : ", ga_instance.generations_completed)
        print("Fitness of the best solution :", ga_instance.best_solution()[1])
        if not (ga_instance.generations_completed % 50):
            #pyRA.FGenome(pyRA.Vect(ga_instance.best_solution()[0])).print()
            pyRA.AGenome(pyRA.Vect(ga_instance.population[ga_instance.best_solution()[2]]), N_REACTIONS).print()

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
                       #parallel_processing=["thread",5],
                       stop_criteria=f"reach_{len(acc_seqs)+len(rej_seqs)}"
                        )

# run the genetic algorithm
ga_problem.run()
np.save("2-parity", best_fitness)
ga_problem.save("2-parity")

# print results
print("Best individual:")
x = pyRA.AGenome(pyRA.Vect(ga_problem.best_solution()[0]), N_REACTIONS)
x.print()
automaton = pyRA.read_agenome(x)

print("Testing on input:")
for seq in acc_seqs[:3]:
    automaton = pyRA.read_agenome(x)
    print("reading:", seq)
    for s in seq:
        print("input: ", s)
        automaton.read_single_input(s)
        automaton.verbose_react()
    while not automaton.is_halted():
        automaton.verbose_react()
    print(f"---------{automaton.accepted()}---------")
for seq in rej_seqs[:3]:
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

# for pure RA
best_fitness = np.zeros(num_generations)
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
                       #parallel_processing=["thread",5],
                       stop_criteria=f"reach_{len(acc_seqs)+len(rej_seqs)}"
                        )

# run the genetic algorithm
ga_problem.run()
np.save("2-parity", best_fitness)
ga_problem.save("2-parity")

# print results
print("Best individual:")
x = pyRA.AGenome(pyRA.Vect(ga_problem.best_solution()[0]), N_REACTIONS)
x.print()
automaton = pyRA.read_agenome(x)

print("Testing on input:")
for seq in acc_seqs[:3]:
    automaton = pyRA.read_agenome(x)
    print("reading:", seq)
    for s in seq:
        print("input: ", s)
        automaton.read_single_input(s)
        automaton.verbose_react()
    while not automaton.is_halted():
        automaton.verbose_react()
    print(f"---------{automaton.accepted()}---------")
for seq in rej_seqs[:3]:
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

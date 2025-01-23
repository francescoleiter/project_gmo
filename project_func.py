
import pyRA

import pygad
import numpy as np


N_REACTIONS = 5

num_generations = 1000
sol_per_pop = 100
num_parents_mating = 70

best_fitness = np.zeros(num_generations)

#           initial state +        5 reactions
num_genes = pyRA.ALPHABET + 3*pyRA.ALPHABET*N_REACTIONS

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

inputs = [pyRA.Vect([a,b,c]) for a in range(0,10) for b in range(0,10) for c in range(0,10)]

def test_func_lin(x):
    return (x[0] + x[1])%3 + 3*x[2]
# UNUSED
def test_func(x):
    return x[0]**2 + x[1]*x[2]

def fitness_function_func(ga_instance, individual, idx, inputs):
    return pyRA.eval_fitness_func(pyRA.FGenome(pyRA.Vect(individual)),test_func, inputs, 40)

def fitness_function_func_pure(ga_instance, individual, idx,  inputs):
    return pyRA.eval_fitness_func_pure(pyRA.FGenome(pyRA.Vect(individual)),test_func, inputs, 40)


def mutation_func_func(offspring, ga_instance):
    for _ in range(int(1/ga_instance.mutation_probability)):
        pyRA.fgenome_mutation_2Darr(offspring, probs, 1.0, p_geom_distribution,
                                 ga_instance.random_mutation_min_val,inh_prob)
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
            pyRA.FGenome(pyRA.Vect(ga_instance.population[ga_instance.best_solution()[2]])).print()

# initialize the class
ga_problem = pygad.GA(
                       gene_type=np.uint16,
                       keep_elitism=keep_elitism,
                       num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       fitness_func=lambda ga_instance, individual, idx: fitness_function_func(ga_instance, individual, idx, inputs=inputs),
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       init_range_low=init_range_low,
                       init_range_high=init_range_high,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       random_mutation_min_val=mutation_min_val,
                       random_mutation_max_val=mutation_max_val,
                       crossover_type=crossover_type,
                       mutation_type=mutation_func_func,
                       mutation_probability=mutation_probability,
                       on_generation=on_gen,
                       #parallel_processing=["thread",5],
                       stop_criteria=f"reach_{0.0}"
                        )

# run the genetic algorithm
ga_problem.run()
np.save("func-lin", best_fitness)
ga_problem.save("func-lin")

# print results
print("Best individual:")
x = pyRA.FGenome(pyRA.Vect(ga_problem.best_solution()[0]))
x.print()
automaton = pyRA.read_fgenome(x)

print("Testing on input:")
for ms in inputs[:6]:
    automaton = pyRA.read_fgenome(x)
    print("reading: ", end = " ")
    ms.print()
    automaton.read_input(ms)
    while not automaton.is_halted():
        automaton.verbose_react()
    print(f"---------{automaton.result_indicator()}---------")

best_fitness = np.zeros(num_generations)
#Using pure automata
ga_problem = pygad.GA(
                       gene_type=np.uint16,
                       keep_elitism=keep_elitism,
                       num_generations=num_generations,
                       num_parents_mating=num_parents_mating,
                       fitness_func=lambda ga_instance, individual, idx: fitness_function_func_pure(ga_instance, individual, idx, inputs=inputs),
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       init_range_low=init_range_low,
                       init_range_high=init_range_high,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       random_mutation_min_val=mutation_min_val,
                       random_mutation_max_val=mutation_max_val,
                       crossover_type=crossover_type,
                       mutation_type=mutation_func_func,
                       mutation_probability=mutation_probability,
                       on_generation=on_gen,
                       #parallel_processing=["thread",5],
                       stop_criteria=f"reach_{0.0}"
                        )

# run the genetic algorithm
ga_problem.run()
np.save("func-lin-pure", best_fitness)
ga_problem.save("func-lin-pure")

# print results
print("Best individual:")
x = pyRA.FGenome(pyRA.Vect(ga_problem.best_solution()[0]))
x.print()

print("Testing on input:")
for ms in inputs[:6]:
    automaton = pyRA.read_fpgenome(x)
    print("reading: ", end = " ")
    ms.print()
    automaton.read_input(ms)
    while not automaton.is_halted():
        automaton.verbose_react()
    print(f"---------{automaton.result_indicator()}---------")
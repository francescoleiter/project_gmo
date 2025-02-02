
import pyRA

import pygad
import numpy as np
# data from https://web.archive.org/web/20060221132749/http://www.nersc.gov/~cding/protein/
data = np.load("data_protein.npy")
print(data.shape)
data_int = list(map(lambda y: list(map(lambda x: int(10 * x), y)), data))
input_data = [pyRA.CVect(line[:-1]) for line in data_int]

N_REACTIONS = 20

num_generations = 500
sol_per_pop = 500
num_parents_mating = 200

best_fitness = np.zeros(num_generations)

#           initial state +        N reactions
num_genes = pyRA.ALPHABET + 3 * pyRA.ALPHABET * N_REACTIONS

# initial values from 0 to 1
init_range_low = 0
init_range_high = 2

# GA parameters
parent_selection_type = "sss"
keep_parents = -1
keep_elitism = 0
crossover_type = None


# improperly used parameters
mutation_max_val = 0.1 # used to set the variance
mutation_probability = 0.9  # used to set the number of mutation to produce
mutation_min_val = 0.1  # probability of being non-zero for reactants and products

inh_prob = 0.02  # probability of an inhibitor to be set on (up to contradictions)
probs = np.array([0.2, 0.8, 0.0])  # probabilities of mutating the initial state, the reactions, the final states


class TestFunc:
    def __init__(self, dataset):
        self.data = dict()
        for line in dataset:
            self.data[tuple(line[:-1])] = line[-1]

    def __call__(self, key):
        return self.data[tuple(key.to_list())]


test_func = TestFunc(data_int)


def fitness_function_func(ga_instance, individual, idx, inputs):
    return pyRA.eval_fitness_func(pyRA.FGenome(pyRA.CVect(individual)), test_func, input_data, 100)


def fitness_function_func_pure(ga_instance, individual, idx, inputs):
    return pyRA.eval_fitness_func_pure(pyRA.FGenome(pyRA.CVect(individual)), test_func, input_data, 100)


def mutation_func_func(offspring, ga_instance):
    pyRA.fgenome_2Dmutation_ngh(offspring, probs, ga_instance.mutation_probability, ga_instance.random_mutation_max_val,
                                ga_instance.random_mutation_min_val, inh_prob)
    return offspring


# function to display the process and update the number of mutations occurring
def on_gen(ga_instance):
    print("variance: ", ga_instance.random_mutation_max_val)
    print("Generation : ", ga_instance.generations_completed)
    print("Fitness of the best solution :", ga_instance.best_solution()[1])
    if not (ga_instance.generations_completed % 10):
        pyRA.FGenome(pyRA.CVect(ga_instance.best_solution()[0])).print()
        # pyRA.AGenome(pyRA.CVect(ga_instance.population[ga_instance.best_solution()[2]]), N_REACTIONS).print()


def on_fit(ga_instance, fitnesses):
    if ga_instance.previous_generation_fitness is None:
        return
    if len(fitnesses[fitnesses <= np.mean(
            ga_instance.previous_generation_fitness[ga_instance.last_generation_parents])]) < 0.2 * len(fitnesses):
        ga_instance.random_mutation_max_val /= 1.1
    else:
        ga_instance.random_mutation_max_val *= 1.1


for index in range(10):

    # initialize the class
    ga_problem = pygad.GA(
                           gene_type=np.uint16,
                           keep_elitism=keep_elitism,
                           num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=lambda ga_instance, individual, idx: fitness_function_func(ga_instance, individual, idx, input_data),
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
                           on_fitness = on_fit,
                           #parallel_processing=["thread",10],
                           stop_criteria=f"reach_{0.0}"

                            )

    # run the genetic algorithm
    ga_problem.run()
    np.save(f"protein_2_{9+index}", best_fitness)
    #ga_problem.save("protein_2")

    # print results
    print("Best individual:")
    x = pyRA.FGenome(pyRA.CVect(ga_problem.best_solution()[0]))
    x.print()


    # for pure RA
    best_fitness = np.zeros(num_generations)
    ga_problem = pygad.GA(
                           gene_type=np.uint16,
                           keep_elitism=keep_elitism,
                           num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=lambda ga_instance, individual, idx: fitness_function_func_pure(ga_instance, individual, idx, input_data),
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
                           on_fitness=on_fit,
                           #parallel_processing=["thread",5],
                           stop_criteria=f"reach_{0.0}"
                            )

    # run the genetic algorithm
    ga_problem.run()
    np.save(f"protein_pure_2_{9+index}", best_fitness)
    #ga_problem.save("protein_pure_4")

    # print results
    print("Best individual:")
    x = pyRA.FGenome(pyRA.CVect(ga_problem.best_solution()[0]))
    x.print()
    automaton = pyRA.read_fpgenome(x)

    print("Testing on input:")

import csv
import random
import copy
import pandas as pd

# this is hyperparameter, need tunning!!
#PPL_SIZE = 50  # the number of DNA in population
#NEW_G_SIZE = 20  # the number od DNA in new generation
#ITERATION = 40  # the number of iteration

def new_nation_n_states(n):
    # this is hyperparameter, need tunning!!
    PPL_SIZE = 50  # the number of DNA in population
    NEW_G_SIZE = 20  # the number od DNA in new generation
    ITERATION = 40  # the number of iteration
    num_ppl_by_states = load_states()
    border_pair = load_borders()
    if n == 1:
        biggest_ppl = num_ppl_by_states.max()
        biggest_state = tuple(num_ppl_by_states[num_ppl_by_states == biggest_ppl].index)
    else:
        ppl = create_starting_population(n, PPL_SIZE, border_pair, num_ppl_by_states)
        for i in range(ITERATION):
            new_g = create_new_generation(ppl, NEW_G_SIZE, border_pair, num_ppl_by_states)
            ppl = pd.concat([ppl,new_g])
            # update new_generation into population
            ppl = ppl.sort_values(by ='ppl', ascending = False)
            ppl = ppl.iloc[:PPL_SIZE,:]
            # survival of the fittest
        biggest_state = ppl.iloc[:1,:].ppl.index[0]
        biggest_ppl = ppl.iloc[:1,:].ppl.values[0]
    return (biggest_state,biggest_ppl)



def load_states():
    states_dict = dict()
    with open('us_states.csv', 'r') as f:
        for line in csv.reader(f):
            states_dict[line[1]] = int(line[3])
    states_series = pd.Series(states_dict)
    states_series.sort_values(ascending = False, inplace = True)

    return states_series


# use dictionary to store the border dict.{'GA': set(['NC','SC'])}. 
def load_borders():
    d = {}
    with open('border_data.csv','r') as f:
        next(f)
        for line in csv.reader(f):
            states = line[1].split('-')
            s1 = states[0]
            s2 = states[1]
            d.setdefault(s1, set()).add(s2)
            d.setdefault(s2, set()).add(s1)
    return d


# create a dataframe of population of different DNA of fixed size N,
# calling function create_new_DNA
# the population would be store in seires
# ('CA','SA')  1000
# dictionary key has to be hashable, which means immutable 
def create_starting_population(n, ppl_size, borders, states_ppl):
    population = {} 
    i = 0
    while i < ppl_size:
        DNA = tuple(create_new_DNA(borders, n)) 
#        if len(DNA) == n and DNA not in population:
        if DNA not in population:
            population[DNA] = fitness(DNA, states_ppl)
            i += 1
            
    return pd.DataFrame.from_dict(population, orient='index', columns=['ppl'])


# DNA is a set, because our states combination is non-ordered
# randomly pick state and then expand/add adjancent state until we have n states in DNA
# make sure the DNA is valid  n-states combination
def create_new_DNA(border_dict, n): 
    DNA = set()
    first_gene = random.choice([*border_dict])
    DNA.add(first_gene) # get a starting state
    gene_candidate_pool = copy.deepcopy(border_dict[first_gene]) # otherwise it will be assigning to the same address
    for i in range(1,n):
        gene = random.choice(list(gene_candidate_pool))
        DNA.add(gene)
        gene_candidate_pool.update(border_dict[gene])
        gene_candidate_pool = gene_candidate_pool - DNA
            
    return sorted(DNA)  # a sorted list

"""
%timeit random.choice([*border_dict])
986 ns ± 19.3 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
%timeit random.choice(list(border_dict))
1.16 µs ± 62.6 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)



%timeit random.choice(list(gene_candidate_pool)) 
878 ns ± 12.4 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
%timeit random.sample(gene_candidate_pool,1)
2.67 µs ± 82.6 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
"""


def fitness(DNA, states_ppl): # tuple DNA
    score = 0
    for gene in DNA:
        score += states_ppl[gene]
    return score

        
# pick parents to give birth to child
# higher fitness, higher chance of being selected
def create_new_generation(population, numofchildren, border_dict, states_ppl):
    tmp = population.ppl
    weights = tmp.values/tmp.sum()
    i = 0
    new_generation = {}
    while i < numofchildren:
        parent = population.sample(n = 1, weights = weights).index[0]  # return tuple
        # this exrtremly slow...1.15ms, if we do equal weights, only 222mu s
        child = mutate(parent, border_dict)  # child is a set
        if validate(child, border_dict):
            child = tuple(sorted(child))
            if child not in new_generation:
                new_generation[child] = fitness(child, states_ppl)
                i += 1
    new_generation = pd.DataFrame.from_dict(new_generation, orient='index', columns=['ppl'])
    return new_generation  #return dictionary of children

    
    
# we only pick DNA to mutate itself valid....
def mutate(DNA, border_dict):  # DNA input as tuple
    # 不需要deepcopy
    valid = False    
    if not valid:
        drop_gene = random.choice(DNA)
        tmp_DNA = set(DNA) - {drop_gene}  # tmp_DNA is a DNA set
        gene_candidate_pool = set()
        for gene in tmp_DNA:
            gene_candidate_pool.update(border_dict[gene])
        gene_candidate_pool = gene_candidate_pool - set(DNA)
        new_gene = random.choice(list(gene_candidate_pool))
        newDNA = tmp_DNA | {new_gene}
    return newDNA  # newDNA output as set


# validate DNA
def validate(DNA, border_dict): # DNA input as set
    valid = True  # assume DNA is valid
    tmp = copy.deepcopy(DNA)
    gene_candidate_pool = set()
    checked_genes = set()
    while valid and len(tmp) != 0:
        gene = tmp.pop()
        gene_candidate_pool.update(border_dict[gene])
        intersection = gene_candidate_pool & tmp
        if len(intersection) == 0: 
            valid = False
        else:
            tmp =  tmp - intersection
            for i in intersection:
                gene_candidate_pool.update(border_dict[gene])
            
    return valid
            

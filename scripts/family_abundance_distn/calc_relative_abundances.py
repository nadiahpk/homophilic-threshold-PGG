import os
import numpy as np
import pandas as pd
from update_functions import add_new, recruit_family

# user parameters
# ---

model_type = 'members_recruit'
model_type = 'leader_driven'
model_type = 'members_attract'

alpha = 3       # weighting to recruit a nonkin
x = 5           # calculate the proportion of the group in the top x most abundant families
mult_n = 2      # when to take recording of abundances, e.g., 10, 100, 1000
n_max = 300000  # maximum group size
q = 0.2         # probability to choose nonkin

res_dir = '../../results/family_abundance_distn/'


# loop through some number of times
# ---

for rep in range(96):

    # collect relative abundances at n, n^2, n^3, ...
    # ---

    # initialise with 1 individual
    n = 1
    abunV = [1] # the abundance classes in the group
    famsV = [1] # the number of families in each abundance class

    rel_abunV = list()
    nV = list()
    while n < n_max:

        n_prev = n
        n = n_prev*mult_n


        # update the abundance and family vectors depending on the model type
        # ---

        if model_type == 'leader_driven':

            for j in range(n_prev+1, n+1): # for each individual

                if np.random.rand() < q:

                    # new individual from new family is recruited

                    if len(abunV) == 1:

                        # if this is the first new individual, append to end
                        # so singletons always in second position
                        abunV.append(1)
                        famsV.append(1)

                    else:

                        # else, add 1 family to the abundance-class 1 (in second position)
                        famsV[1] += 1

                else:

                    # new individual is recruited to leader's family (assumed the first one)
                    abunV[0] += 1

        elif model_type == 'members_recruit':

            for j in range(n_prev+1, n+1): # for each individual

                if np.random.rand() < q:

                    # new individual from new family is recruited
                    abunV, famsV = add_new(abunV, famsV)

                else:

                    # new individual is recruited by a randomly chosen current group member
                    abunV, famsV = recruit_family(abunV, famsV)

        else: # model_type == 'members_attract':

            for j in range(n_prev+1, n+1): # for each individual

                if np.random.rand() < alpha / (alpha+j-1):

                    # new individual from new family is recruited
                    abunV, famsV = add_new(abunV, famsV)

                else:

                    # new individual is recruited by a randomly chosen current group member
                    abunV, famsV = recruit_family(abunV, famsV)

        # calculate the proportion of individuals in the x most abundant families
        # ---

        abunT, famsT  = zip(*sorted(zip(abunV, famsV), key=lambda v: v[0], reverse=True))
        abunV = list(abunT)
        famsV = list(famsT)

        # create a list of the abundances in the top x families

        top_abunV = list()
        cnt = 0
        for abun, fams in zip(abunV, famsV):
            
            top_abunV += [abun]*fams
            cnt += fams

            if cnt >= x:
                break

        # fix the length of top_abunV
        len_top_abunV = len(top_abunV)
        if len(top_abunV) > x:
            top_abunV = top_abunV[:x]
        elif len(top_abunV) < x:
            top_abunV += [0]*(x-len_top_abunV)

        rel_abun = [ top_abunV[i] / n for i in range(x) ]

        # store
        # ---

        rel_abunV.append(rel_abun)
        nV.append(n)


    # for each abundandance class, write a row of results to csv
    # ---

    columns = ['n=' + str(n) for n in nV ]
    for xi in range(x):

        rel_abun_x = [ rel_abun[xi] for rel_abun in rel_abunV ]
        df = pd.DataFrame.from_records([rel_abun_x], columns=columns)

        fname = res_dir + 'relative_abundances_' + model_type + '_rank_' + str(xi+1) + '.csv'

        if not os.path.isfile(fname):
            df.to_csv(fname, mode='w', header=True, index=False) # add with header
        else:
            df.to_csv(fname, mode='a', header=False, index=False) # append

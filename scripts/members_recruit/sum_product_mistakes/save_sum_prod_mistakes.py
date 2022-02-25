# Save the product of mistake indices, which are used to calculate the
# partition probabilities when there is a fixed probability of 
# mistakenly choosing a stranger instead of a family member, q.

import pandas as pd
import itertools as it
from sympy.utilities.iterables import multiset_permutations
from scipy.special import comb
import numpy as np

import sys
sys.path.append('../../../functions/')

from my_functions import read_matrix_M


# parameters
# ---

group_sizes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18] 

# The limit is n = 18. Size 19 returns error:
#   doing 22 of 490
#   RuntimeWarning: overflow encountered in long_scalars
#   sum_prod_mistakes += count_orderings*prod_mistake


# for each n, find the product of mistake probabilities and save
# ---

for n in group_sizes:

    print('-----------')
    print('doing n = ' + str(n))


    # partns lists all the possible partitions of n
    # ---

    lm, partns, M_num, M_den = read_matrix_M('../../../results/matrix_M/matrix_M' + str(n) + '.csv')
    num_partns = len(partns)


    # create strings of all possible partitions (use later as rows in the .csv table)
    # ---

    partn_strings = [ '|'.join([ str(partn_ij) for partn_ij in partn_i]) for partn_i in partns ]


    # for each partition in partns, calculate sum_prod_mistakes
    # ---

    sum_prod_mistakesV = list()
    for idx, partn in enumerate(partns):

        print('doing partition ' + str(idx) + ' of ' + str(num_partns))

        # for each possible ordering of families joining the group
        # find the sum product term involving the mistake indices

        Phi = len(partn)           # the total number of families
        sum_prod_mistakes = 0   # initialise sum product term
        for nV in multiset_permutations(partn):

            # create a list of possible mistake indices vectors and their products

            max_ms = np.cumsum(nV[:-1]) # find maximum mistake index
            mVs = list(filter( lambda mV: all( mi <= m_max for mi, m_max in zip(mV, max_ms) ), it.combinations(range(1, n), Phi-1) ))
            prod_mistakes = [ np.prod(mV) for mV in mVs ]

            # for each possible mistake indices vector, 
            # given this ordering of families joining the group,
            # calculate how many orderings of individuals joining the group correspond to that

            for mV, prod_mistake in zip(mVs, prod_mistakes):

                # calculate how many orderings correspond
                count_orderings = 1
                for j in range(2, Phi+1):

                    num = n - mV[j-2] - 1 - sum(nV[j:Phi+1])    # n - m_{j-1} - 1 - sum_{k=j+1}^Phi n_{i_k}
                    den = nV[j-1] - 1                           # n_{i_j} - 1
                    count_orderings *= comb(num, den, exact=True)

                # so there are count_orderings arrival orderings that correspond to
                # this mistake index vector and its product
                sum_prod_mistakes += count_orderings*prod_mistake

        sum_prod_mistakesV.append(int(sum_prod_mistakes))

        # this is how I would calculate the probability of this outcome:
        #   P = (np.prod([factorial(ni-1) for ni in partn]) / factorial(n-1)) * q**(Phi-1) * (1-q)**(n-Phi) * sum_prod_mistakes


    # save it to a csv file
    # ---

    # two columns, the partition, and the product of mistake indices term
    df = pd.DataFrame(list(zip(partn_strings, sum_prod_mistakesV)), 
            columns =['partition', 'sum_product_mistake_indices'])

    # write
    df.to_csv('../../../results/members_recruit/sum_product_mistakes/sum_prod_mistakes' + str(n) + '.csv', index=False)


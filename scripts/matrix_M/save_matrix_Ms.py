# For convenience, save in csv files the numerator and denominator of
# the matrices M that convert from the probabilities of possible
# partitions F_{n -> lambda} to the n-relatednesses theta_{l -> m}
# e.g., the matrix in Eq. C1 in Ohtsuki (2014)

import itertools as it
from sympy.utilities.iterables import multiset_permutations
from scipy.special import comb
import numpy as np
import pandas as pd
import os

import sys
sys.path.append("../../functions/")

from my_functions import partitionInteger, read_matrix_M


# parameters
# ---

group_sizeV = [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 ]

for n in group_sizeV:

    print('-----------------------------------')
    print('DOING n = ' + str(n))


    # get old answer to check correctness of this code
    # ---

    # this is where I stored answers from an earlier, slower algorithm
    fname = '../../../2109_genetic/bk_eql_0/results/matrix_M' + str(n) + '.csv'

    if os.path.isfile(fname):

        can_check = True
        lm, F, M_num_orig, M_den_orig = read_matrix_M(fname)

    else:

        can_check = False


    # get list of (l,m) values for theta_{l->m} and lists of partitions needed
    # ---

    # create list of (l,m), the possible theta_{l -> m} indices, i.e., the lhs of Eq C1
    lm = [ (l, m) for l in range(1, n+1) for m in range(1, l+1) ]
    # e.g., [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (4, 3), (4, 4)]

    # get the partitions of every integer 2, .., n
    # we'll use later to create our list of q_i
    P = [ list(partitionInteger(v)) for v in range(2, n+1) ]

    # get every partition of the group size n, i.e., the column vector on rhs of Eq C1
    F = P[-1]
    # e.g., F = [[1, 1, 1, 1], [1, 1, 2], [1, 3], [2, 2], [4]]

    print('found partitions')


    # create the matrix M
    # ---

    # size of matrix M
    no_rows = len(lm)
    no_cols = len(F)


    # create the denominator matrix for M
    # ---

    # the values of the denominator depend on l only
    M_denD = { l: comb(n, l, exact=True) for l in range(1, n+1) }

    # the first row corresponds to theta_{1->1}, which is all ones
    M_den = np.zeros((no_rows, 1), dtype=int)
    M_den[0, :] = 1

    # place the denominator value in each corresponding row

    for row_, (l, m) in enumerate(lm[1:]):

        M_den[row_ + 1, :] = M_denD[l]


    # create the numerator matrix for M
    # ---

    # default entry is 0
    M_num = np.zeros((no_rows, no_cols), dtype=int)

    # the first row corresponds to theta_{1->1}, which is all ones
    M_num[0, :] = 1

    # the ending rows correspond to theta_{n->m}, which is 1 wherever m = w
    wV = [ len(nV) for nV in F ]
    for row in range(no_rows - n, no_rows):

        m = lm[row][1]
        cols = [ i for i, w in enumerate(wV) if w == m ]

        # put a 1 wherever the condition m = w satisfied
        M_num[row, cols] = 1

    # for the rest, we must calculate ...


    # loop through each row of M's numerator and calculate entries' values
    # ---

    for row_, (l, m) in enumerate(lm[1:-n]):

        print('doing (l,m) = ' + str(l) + ', ' + str(m) )


        # find qVs, which are the possible (q1, q2,..,qm), where where q_i is the number of individuals in family i
        # ---

        # pull out partitions of l with length m (indexing starts at 2)
        P_lm = [ PV for PV in P[l-2] if len(PV) == m ]

        # find every multiset permutation of each partition in P_lm
        # e.g., [1, 1, 1, 4] -> [[1, 1, 1, 4], [1, 1, 4, 1], [1, 4, 1, 1], [4, 1, 1, 1]]

        qVsV = list()

        for PV in P_lm:

            qVsV.append(list(multiset_permutations(PV)))

        
        # for each column, calculate the number of ways of choosing l individuals from m families
        # ---

        for col, (w, nV) in enumerate(zip(wV,F)):

            # nV is the partition, i.e., the entry in the column vector on the rhs of Eq C1
            # w is the number of families in this partition i.e., the partition's length

            if m <= w: # we can't choose m families if there aren't enough in this partition

                # to save time, precheck which PV in P_lm is even possible with this nV
                possibleV = [ all([ ni >= Pi for ni, Pi in zip(reversed(nV), reversed(PV)) ]) for PV in P_lm ]

                if any(possibleV):

                    # find iVs, which is a list of indices of which families we draw from
                    iVs = list(it.combinations(range(w), m))

                    # take the sum over i and over q_i (main equation)
                    M_num_v = 0

                    # this is to skip multiset permutations of partitions of l of length m that are not possible anyway
                    for qVs, possible in zip(qVsV, possibleV):
                        if possible:

                            # only do the expensive loop for those that are possible
                            for iV in iVs:
                                for qV in qVs:

                                    if all([ qi <= nV[i] for qi, i in zip(qV, iV) ]):

                                        M_num_v += np.prod([ comb(nV[i], qi, exact=True) for i, qi in zip(iV,qV) ])

                    # store result in matrix
                    M_num[row_ + 1, col] = M_num_v


    # check it matches the previous result (if we've got one)
    # ---

    if can_check:

        check_bool = (M_num_orig == M_num).all()

        if check_bool:

            print('New algorithm result matches previous')

        else:

            print('New algorithm result DOES NOT MATCH previous')


    # write to csv file
    # ---

    # create the first two columns, which are the theta indices
    df_lm = pd.DataFrame(lm, columns = ['theta_idx_l', 'theta_idx_m'])

    # the numerator of M
    df_num = pd.DataFrame(M_num, columns = [ 'num_' + '|'.join([ str(Fij) for Fij in Fi]) for Fi in F ])
    # reverse column headers with: ss.split('_')[1].split('|')

    # the denominator of M
    df_den = pd.DataFrame(M_den, columns = ['den'])

    # merge them
    df = pd.concat([df_lm, df_den, df_num], axis=1)

    # write
    df.to_csv('../../results/matrix_M/matrix_M' + str(n) + '.csv', index=False)

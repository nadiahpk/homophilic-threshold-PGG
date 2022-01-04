from scipy.special import binom
import itertools as it
from math import factorial
import numpy as np
import pandas as pd

def partitionInteger(n):
    '''
    Copied verbatim from: http://jeromekelleher.net/generating-integer-partitions.html
    Finds all partitions of the integer n
    '''

    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def read_matrix_M(fname):
    '''
    Read the saved numerator and denominator of the matrix M (converts 
    partition probabilities F_{n -> lambda} to the n-relatednesses 
    theta_{l -> m}). Assumes the .csv file was created by save_matrix_Ms.py

    Inputs
    ---

    fname, str
        The name of the csv file where the matrix is saved

    Outputs
    ---

    lm, list of tuples of ints
        A list of (l,m), the possible theta_{l -> m} indices

    F, list of lists of ints
        Every partition of n

    M_num, np matrix of ints
        The numerator of matrix M

    M_den, np matrix of ints
        The denominator of matrix M
    '''

    # read in the dataframe where it is stored
    df = pd.read_csv(fname)

    # the first two columns give (l,m), the indices of theta_{l -> m}
    lm = list(zip(df['theta_idx_l'].values, df['theta_idx_m']))

    # the headers give the F partitions
    F_str = [ s.split('_')[1] for s in df.columns if s[:3] == 'num'] # use the numerator one
    F = [ [ int(si) for si in s.split('|')] for s in F_str ]

    # recreate the M_num matrix
    M_num = np.array(df[ ['num_' + s for s in F_str] ])

    # recreate the M_den matrix

    # check if we have stored the entire denominator or just the first column (for compatibility with old code)
    if 'den' in df.columns:     # stored first column only

        M_den_col = np.matrix(df['den'])
        M_den = np.array(np.tile(M_den_col.transpose(), (1, len(F_str)))) # make array again for backward compatibility

    else:                       # stored entire matrix

        M_den = np.array(df[ ['den_' + s for s in F_str] ])

    return(lm, F, M_num, M_den)


# ---------------------------------------------------------------------
# TODO - document below
def get_PV_leader(n, F, q):

    P = list()
    for Fi in F:

        if len(Fi) == 1:

            Pi = (1-q)**(n-1)

        else:

            # rely on the largest group being at the end in code below
            Fi.sort()

            if Fi[-2] != 1:

                # this cannot happen in the leader model because all strangers are unrelaeted
                Pi = 0

            else:

                s = len(Fi)-1
                Pi = binom(n-1, s) * q**s * (1-q)**(n-1-s)

        P.append(Pi)

    return(P)

def get_PV_constant_mistake(n, F, q, fname): 
    '''

    Inputs
    ---

    n, int
        Group size

    F, list of list of ints
        A list of partitions of n

    q, float
        Probability to mistakenly choose a non-family member

    fname, string
        Location of the sum_prod_mistakes[n].csv file

    '''

    # deal with special cases
    # ---

    if q == 0:

        # never mistakenly choose a stranger
        PV = [ 1 if Fi == [n] else 0 for Fi in F ]

    elif q == 1: 

        # always choose a stranger
        PV = [ 1 if Fi == [1]*n else 0 for Fi in F ]

    else:

        # read in and prepare the needed info
        # ---

        df = pd.read_csv(fname)
        df.set_index('partition', inplace=True)


        # for each partition, calculate the probability
        # ---

        PV = list()
        for Fi in F:

            # get sum_prod_mistakes and other info

            Phi = len(Fi) # the total number of families
            Fi_str = '|'.join([ str(Fij) for Fij in Fi]) 
            sum_prod_mistakes = df.loc[Fi_str]['sum_product_mistake_indices']

            # make calculation

            P = (np.prod([factorial(phi-1) for phi in Fi]) / factorial(n-1)) * q**(Phi-1) * (1-q)**(n-Phi) * sum_prod_mistakes
            PV.append(P)

    return(PV)


def get_PV_stranger_weighting(n, F, alpha):

    # place to store probabilities of each partition
    P = [0]*len(F)

    fact_n = factorial(n)
    for i, partn in enumerate(F):

        Phi = len(partn)

        # phi_i i s the number of species with i individuals
        # create a list of non-zero (i, phi_i) pairs

        iV = set(partn) # which i's occur
        i_phiV = [ (i, partn.count(i)) for i in iV ]

        # get terms in denom
        AA = np.prod([ i**phi_i for i, phi_i in i_phiV ])
        BB = np.prod([ factorial(phi_i) for _, phi_i in i_phiV ])
        CC = np.prod([ alpha + k - 1 for k in range(1, n+1) ])

        # calc P
        P[i] = fact_n * alpha**Phi / (AA*BB*CC)

    return(P)



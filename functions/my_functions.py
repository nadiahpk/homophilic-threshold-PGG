from scipy.special import binom
import itertools as it
from math import factorial
import numpy as np
import pandas as pd


def partitionInteger(n):
    '''
    Finds all partitions of the integer n.

    Copied verbatim from: http://jeromekelleher.net/generating-integer-partitions.html

    Kelleher, J. and O’Sullivan, B. (2009). Generating all partitions: 
    a comparison of two encodings, arXiv preprint arXiv:0909.2331
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

    partns, list of lists of ints
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

    # the headers give the partitions
    partn_strs = [ s.split('_')[1] for s in df.columns if s[:3] == 'num'] # use the numerator one
    partns = [ [ int(si) for si in s.split('|')] for s in partn_strs ]

    # recreate the M_num matrix
    M_num = np.array(df[ ['num_' + s for s in partn_strs] ])

    # recreate the M_den matrix

    # check if we have stored the entire denominator or just the first column 
    # (for compatibility with old code)
    if 'den' in df.columns:     # stored first column only

        M_den_col = np.matrix(df['den'])

        # make array again for backward compatibility
        M_den = np.array(np.tile(M_den_col.transpose(), (1, len(partn_strs))))

    else:                       # stored entire matrix

        M_den = np.array(df[ ['den_' + s for s in partn_strs] ])

    return(lm, partns, M_num, M_den)


# functions for random group formation
# ---

def calc_p_peak_random_group(parD):
    '''
    Find the value of p that maximises Delta p when groups are formed 
    randomly (no homophily)

    Inputs
    ---

    parD, dict
        Parameter values with keys below
        - 'tau': int, Min no. Cooperators for public good productn
        - 'n': int, Group size
        - 'W': float, Cooperator payoff if threshold met
        - 'X': float, Cooperator payoff if threshold not met
        - 'Y': float, Defector payoff if threshold met
        - 'Z': float, Defector payoff if threshold not met


    Outputs
    ---

    p_peak, float [0,1]
        Value of p that maximises Delta p
    '''

    # unpack dictionary of parameter values
    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']

    # find the p value that maximises Delta p
    p_peak = (tau-1)*(W-X)/((tau-1)*(W-X)+(Y-Z)*(n-tau))

    return p_peak


# functions for homophilic group formation
# ---

# obtain F from degree of homophily

def get_FV_leader_probability(n, partns, q):
    '''
    Using the 'leader recruits + nonkin probability' homophilic 
    group-formation model, calculate the vector of partition 
    probabilities F.


    Inputs:
    ---

    n, int
        Group size

    partns, list of list of ints
        List of all possible partitions of n

    q, float (0,1)
        The probability that the recruitor makes a mistake and 
        recruits a nonkin new member


    Outputs:
    ---

    FV, list of floats
        The probability that a group is formed with family partition 
        structure corresponding to structures in partns
    '''

    FV = list()
    for partn in partns:

        if len(partn) == 1:

            Fi = (1-q)**(n-1)

        else:

            # rely on the largest group being at the end in code below
            partn.sort()

            if partn[-2] != 1:

                # this cannot happen in the leader model because all strangers are unrelaeted
                Fi = 0

            else:

                s = len(partn)-1
                Fi = binom(n-1, s) * q**s * (1-q)**(n-1-s)

        FV.append(Fi)

    return FV 


def get_FV_random_probability(n, partns, q, fname=None): 
    '''
    Using the 'random recruitor + nonkin probability' homophilic 
    group-formation model, calculate the vector of partition 
    probabilities F.

    Inputs:
    ---

    n, int
        Group size

    partns, list of list of ints
        List of all possible partitions of n

    q, float (0,1)
        The probability that the recruitor makes a mistake and 
        recruits a nonkin new member

    fname, string
        Location of the sum_prod_mistakes[n].csv file


    Outputs:
    ---

    FV, list of floats
        The probability that a group is formed with family partition 
        structure corresponding to structures in partns
    '''

    if fname is None:
        # where sum_prod_mistakes stored
        fname = '../../results/random_recruitor_nonkin_probability/sum_product_mistakes/sum_prod_mistakes' + str(n) + '.csv' 

    # deal with special cases
    # ---

    if q == 0:

        # never mistakenly choose a stranger
        FV = [ 1 if partn == [n] else 0 for partn in partns ]

    elif q == 1: 

        # always choose a stranger
        FV = [ 1 if partn == [1]*n else 0 for partn in partns ]

    else:

        # read in and prepare the needed info
        # ---

        df = pd.read_csv(fname)
        df.set_index('partition', inplace=True)


        # for each partition, calculate the probability
        # ---

        FV = list()
        for partn in partns:

            # get sum_prod_mistakes and other info

            Phi = len(partn) # the total number of families
            partn_str = '|'.join([ str(partn_ij) for partn_ij in partn]) 
            sum_prod_mistakes = df.loc[partn_str]['sum_product_mistake_indices']

            # make calculation

            F = (np.prod([factorial(phi-1) for phi in partn]) / factorial(n-1)) * q**(Phi-1) * (1-q)**(n-Phi) * sum_prod_mistakes
            FV.append(F)

    return FV


def get_FV_random_weighting(n, partns, alpha):
    '''
    Using the 'random recruitor + nonkin weighting' homophilic 
    group-formation model, calculate the vector of partition 
    probabilities F.


    Inputs:
    ---

    n, int
        Group size

    partns, list of list of ints
        List of all possible partitions of n

    alpha, float (0,infty)
        The 'stranger weighting', the weighting given to recruitment 
        of a nonkin new member.


    Outputs:
    ---

    FV, list of floats
        The probability that a group is formed with family partition 
        structure corresponding to structures in partns
    '''


    # place to store probabilities of each partition
    FV = [0]*len(partns)

    fact_n = factorial(n)
    for i, partn in enumerate(partns):

        Phi = len(partn)

        # phi_i i s the number of species with i individuals
        # create a list of non-zero (i, phi_i) pairs

        iV = set(partn) # which i's occur
        i_phiV = [ (i, partn.count(i)) for i in iV ]

        # get terms in denom
        AA = np.prod([ i**phi_i for i, phi_i in i_phiV ])
        BB = np.prod([ factorial(phi_i) for _, phi_i in i_phiV ])
        CC = np.prod([ alpha + k - 1 for k in range(1, n+1) ])

        # calc F
        FV[i] = fact_n * alpha**Phi / (AA*BB*CC)

    return FV

def calc_delta_p_random_probability(parD, lm, partns, M, q, p):
    '''
    Using the 'random recruitor + nonkin probability' homophilic 
    group-formation model, calculate delta_p, which is proportional to 
    the change in frequency of Cooperators in the population.

    Inputs
    ---

    parD, dict
        Parameter values with keys below
        - 'tau': int, Min no. Cooperators for public good productn
        - 'n': int, Group size
        - 'W': float, Cooperator payoff if threshold met
        - 'X': float, Cooperator payoff if threshold not met
        - 'Y': float, Defector payoff if threshold met
        - 'Z': float, Defector payoff if threshold not met

    lm, list of tuples of ints
        A list of (l,m), the possible theta_{l -> m} indices

    partns, list of list of ints
        List of all possible partitions of n

    M, matrix of floats
        Matrix that converts partition probabilities F_{n -> lambda} 
        to the n-relatednesses theta_{l -> m}). 

    q, float (0,1)
        The probability that the recruitor makes a mistake and 
        recruits a nonkin new member
        
    p, float (0,1)
        The proportion of Cooperators in the population


    Outputs
    ---

    delta_p, float
        A value that is proportional to the change in frequency of 
        Cooperators in the population.
    '''

    # unpack parameter values
    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']

    # payoff functions
    b = lambda k: Z if k < tau else Y
    a = lambda k: X if k < tau-1 else W

    # find the probability of each partition according to the stranger-weighting model
    FV = get_FV_random_probability(n, partns, q)

    # calculate theta_{l -> m}
    thetaV = M @ np.array(FV)

    # create a dictionary from (l,m) to theta_{l -> m}
    thetaD = dict(zip(lm, thetaV))

    # create the rho function
    rho = lambda l: 1 if l == 0 else sum(thetaD[(l,m)]*p**m for m in range(1, l+1))

    # calculate Delta p
    delta_p = sum( 
            sum( (-1)**(l-k) * binom(l, k) * binom(n-1, l) * 
                ( (1-rho(1))*rho(l+1)*a(k) - rho(1)*(rho(l)-rho(l+1))*b(k) ) 
                for l in range(k, n) ) 
            for k in range(0, n) )

    return delta_p

def q0_root_eqn_random_probability(parD, lm, partns, M, q):
    '''
    Using the 'random recruitor + nonkin probability' homophilic 
    group-formation model.

    Solve this equation to calculate q_0, the mistake probability at 
    which the unstable isocline meets the p=0 axis.

    Inputs
    ---

    parD, dict
        Parameter values with keys below
        - 'tau': int, Min no. Cooperators for public good productn
        - 'n': int, Group size
        - 'W': float, Cooperator payoff if threshold met
        - 'X': float, Cooperator payoff if threshold not met
        - 'Y': float, Defector payoff if threshold met
        - 'Z': float, Defector payoff if threshold not met

    partns, list of list of ints
        List of all possible partitions of n

    M, matrix of floats
        Matrix that converts partition probabilities F_{n -> lambda} 
        to the n-relatednesses theta_{l -> m}). 

    q, float (0,1)
        The probability that the recruitor makes a mistake and 
        recruits a nonkin new member
        

    Outputs
    ---

    res, float
        Solve res=0 to obtain q_0
    '''

    # unpack parameter values
    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']

    # payoff functions
    b = lambda k: Z if k < tau else Y
    a = lambda k: X if k < tau-1 else W
    
    # find the probability of each partition according to the stranger-weighting model
    FV = get_FV_random_probability(n, partns, q)

    # calculate theta_{l -> m}
    thetaV = M @ np.array(FV)

    # create a dictionary from (l,m) to theta_{l -> m}
    thetaD = dict(zip(lm, thetaV))
    res = X + (W-X)*sum( binom(n-1, l) * thetaD[(l+1, 1)] * (-1)**(l-tau+1) * binom(l-1,tau-2) 
            for l in range(tau-1, n) ) - b(0)

    return res

# find the root of this equation to find where the isocline intercepts with the p_0 = 1 axis
def q1_root_eqn_random_probability(parD, lm, partns, M, q):
    '''
    Using the 'random recruitor + nonkin probability' homophilic 
    group-formation model.

    Solve this equation to calculate q_1, the mistake probability at 
    which the stable isocline meets the p=1 axis.

    Inputs
    ---

    parD, dict
        Parameter values with keys below
        - 'tau': int, Min no. Cooperators for public good productn
        - 'n': int, Group size
        - 'W': float, Cooperator payoff if threshold met
        - 'X': float, Cooperator payoff if threshold not met
        - 'Y': float, Defector payoff if threshold met
        - 'Z': float, Defector payoff if threshold not met

    partns, list of list of ints
        List of all possible partitions of n

    M, matrix of floats
        Matrix that converts partition probabilities F_{n -> lambda} 
        to the n-relatednesses theta_{l -> m}). 

    q, float (0,1)
        The probability that the recruitor makes a mistake and 
        recruits a nonkin new member
        

    Outputs
    ---

    res, float
        Solve res=0 to obtain q_1
    '''

    # unpack parameter values
    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']

    # payoff functions
    b = lambda k: Z if k < tau else Y
    a = lambda k: X if k < tau-1 else W
    
    # find the probability of each partition according to the stranger-weighting model
    FV = get_FV_random_probability(n, partns, q)

    # calculate theta_{l -> m}
    thetaV = M @ np.array(FV)

    # create a dictionary from (l,m) to theta_{l -> m}
    thetaD = dict(zip(lm, thetaV))
    res = Y + (Z-Y)*sum( binom(n-1, l) * thetaD[(l+1, 1)] * (-1)**(l-(n-tau)) * binom(l-1,n-tau-1) 
            for l in range(n-tau, n) ) - W

    return res

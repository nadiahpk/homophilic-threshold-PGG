from math import factorial
from scipy.special import binom
import numpy as np
from scipy.optimize import bisect
import pandas as pd

import sys
sys.path.append('../../functions/')

from my_functions import read_matrix_M, calc_p_peak_random_group, calc_delta_p_members_attract, alpha0_root_eqn_members_attract, alpha1_root_eqn_members_attract

# ----------------------

# user parameters
# ---

res_dir = '../../results/members_attract/'

# how many increments on the x-axis for plotting
ngrid = 50

nV =   [4, 6, 8, 10, 12, 14, 16, 18]
tauV = [2, 3, 4,  5,  6,  7,  8,  9]

# which default parameter values to use from the file default_params.csv
suffix = '_6'
suffix = '_7'
suffix = '_5'


# fixed parameters
# ---

# all parameter names corresponding to columns in default_params.csv
par_names = ['n', 'tau', 'W', 'X', 'Y', 'Z']


# read in the default parameter values
# ---

df = pd.read_csv('../default_params.csv')
df.set_index('suffix', inplace=True)
parD = { par_name: df[par_name].loc[suffix] for par_name in par_names }

# unpack parameter values
W = parD['W']
X = parD['X']
Y = parD['Y']
Z = parD['Z']


# find how the critical values vary with n
# ---

# places to store results
alpha0V = list()
alpha1V = list()
alphahV = list()   # store alpha_star values

for n, tau in zip(nV, tauV):

    print('doing n = ' + str(n) + ', tau = ' + str(tau))

    # update parameter values
    parD['n'] = n
    parD['tau'] = tau

    # get previously saved lm, partns, and M
    # ---

    lm, partns, M_num, M_den = read_matrix_M('../../results/matrix_M/matrix_M' + str(n) + '.csv')
    M = M_num / M_den # the matrix for converting F to theta


    # find alpha_0, the maximum mistake probability where collectivists can still invade
    # ---

    # should always be able to invade when perfect homophily
    alpha_lo = 0
    # assert np.sign(alpha0_root_eqn_members_attract(parD, lm, partns, M, alpha_lo)) == 1

    # should never be able to invade when randomly formed
    # assert np.sign(alpha0_root_eqn_members_attract(parD, lm, partns, M, np.inf)) == -1

    # find an upper bound on alpha where cooperators can no longer invade by doubling
    alpha_hi = 1
    while np.sign(alpha0_root_eqn_members_attract(parD, lm, partns, M, alpha_hi)) == 1: # while C can invade
        alpha_hi *= 2

    # now you can find the alpha_0 between them
    alpha_0 = bisect(lambda alpha: alpha0_root_eqn_members_attract(parD, lm, partns, M, alpha), alpha_lo, alpha_hi) # NOTE
    alpha0V.append(alpha_0)


    # find alpha_1, the minimum mistake probability where individualists can invade
    # ---

    # first, verify a high and low
    alpha_lo = 0

    # find an upper bound by doubling
    alpha_hi = 1
    while np.sign(alpha1_root_eqn_members_attract(parD, lm, partns, M, alpha_hi)) == -1: # while D cannot invade
        alpha_hi *= 2

    # find alpha_1 where D cannot invade
    alpha_1 = bisect(lambda alpha: alpha1_root_eqn_members_attract(parD, lm, partns, M, alpha), alpha_lo, alpha_hi)
    alpha1V.append(alpha_1)


    # find alpha_hat, if it exists
    # ---

    # find p_peak (provides bound for searching for isoclines)
    p_peak = calc_p_peak_random_group(parD)


    # find the maximum alpha where there are two isoclines

    # if Delta p at p_peak < 0, then there are zero interior equilibria when alpha=1;
    # otherwise, there are two
    delta_p_at_alphainf = calc_delta_p_members_attract(parD, lm, partns, M, np.inf, p_peak)

    if delta_p_at_alphainf > 0:

        alpha_max = np.inf
        p_sing = p_peak

    else:

        # find where the transition from 2 to 0 equilibria is
        # ---

        # first, find an alpha_hi where we lose the two interior equilibria 

        alpha_hi = max((alpha_0, alpha_1))*2
        while calc_delta_p_members_attract(parD, lm, partns, M, alpha_hi, p_peak) > 0:
            alpha_hi *= 2

        # now find the point between alpha_lo and alpha_hi where the transition happens

        alpha_lo = alpha_1 if alpha_1 > alpha_0 else alpha_0
        alpha_mid = (alpha_lo + alpha_hi)/2
        p_mid = p_peak

        while abs(alpha_lo-alpha_hi) > 1e-4:

            # try to find the upper and lower isoclines at alpha_mid

            try:
                p_up = bisect(lambda p: calc_delta_p_members_attract(parD, lm, partns, M, alpha_mid, p), 
                        p_mid, 1-1e-6)
            except:
                p_up = np.nan

            try:
                p_lo = bisect(lambda p: calc_delta_p_members_attract(parD, lm, partns, M, alpha_mid, p), 
                        1e-6, p_mid)
            except:
                p_lo = np.nan

            if np.isnan(p_lo) or np.isnan(p_lo):
                alpha_hi = alpha_mid                # search to the left
                alpha_mid = (alpha_lo + alpha_mid)/2
            else:
                alpha_lo = alpha_mid                # search to the right
                alpha_mid = (alpha_hi + alpha_mid)/2
                p_mid = (p_up + p_lo)/2

        alpha_max = alpha_lo
        p_sing = p_mid

    # store
    alphahV.append(alpha_max)


# write isocline data to file
# ---


# nV, alpha0V, alpha1V, alphahV
df_out = pd.DataFrame(list(zip(nV, alpha0V, alpha1V, alphahV)), columns=['nV', 'alpha0V', 'alpha1V', 'alphahV'])
fname = res_dir + 'critical_v_n' + suffix + '.csv'
df_out.to_csv(fname, index=False)

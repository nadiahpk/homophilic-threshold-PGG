# create a detailed example

import os
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

# how many increments on the x-axis for plotting
ngrid = 50

# which default parameter values to use from the file default_params.csv
suffix = '_6'
suffix = '_5'
suffix = '_7'

res_dir = '../../results/members_attract/'

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
n = parD['n']
tau = parD['tau']
W = parD['W']
X = parD['X']
Y = parD['Y']
Z = parD['Z']


# get previously saved lm, partns, and M
# ---

lm, partns, M_num, M_den = read_matrix_M('../../results/matrix_M/matrix_M' + str(n) + '.csv')
M = M_num / M_den # the matrix for converting F to theta


# find alpha_0, the maximum mistake weighting where collectivists can still invade
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
alpha_0 = bisect(lambda alpha: alpha0_root_eqn_members_attract(parD, lm, partns, M, alpha), alpha_lo, alpha_hi)
print('alpha_0 = ' + str(alpha_0))


# find alpha_1, the minimum mistake weighting where individualists can invade
# ---

# first, verify a high and low
alpha_lo = 0

# find an upper bound by doubling
alpha_hi = 1
while np.sign(alpha1_root_eqn_members_attract(parD, lm, partns, M, alpha_hi)) == -1: # while D cannot invade
    alpha_hi *= 2

# find alpha_1 where D cannot invade
alpha_1 = bisect(lambda alpha: alpha1_root_eqn_members_attract(parD, lm, partns, M, alpha), alpha_lo, alpha_hi)
print('alpha_1 = ' + str(alpha_1))

# find p_peak (provides bound for searching for isoclines)
# ---

p_peak = calc_p_peak_random_group(parD)


# find the maximum alpha where there are two isoclines
# ---

# if Delta p at p_peak is < 0, then there are zero interior equilibria when alpha=inf;
# otherwise, there are two
delta_p_at_alphainf = calc_delta_p_members_attract(parD, lm, partns, M, np.inf, p_peak)

if delta_p_at_alphainf > 0:

    alpha_hat = np.inf
    p_at_alpha_hat = p_peak

else:

    # find where the transition from 2 to 0 equilibria is
    # ---

    print('-----')
    print('finding where the transition from 2 to 0 equilibria is')

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
            p_up = bisect(lambda p: calc_delta_p_members_attract(parD, lm, partns, M, alpha_mid, p), p_mid, 1-1e-6)
        except:
            p_up = np.nan

        try:
            p_lo = bisect(lambda p: calc_delta_p_members_attract(parD, lm, partns, M, alpha_mid, p), 1e-6, p_mid)
        except:
            p_lo = np.nan

        if np.isnan(p_lo) or np.isnan(p_lo):
            alpha_hi = alpha_mid                # search to the left
            alpha_mid = (alpha_lo + alpha_mid)/2
        else:
            alpha_lo = alpha_mid                # search to the right
            alpha_mid = (alpha_hi + alpha_mid)/2
            p_mid = (p_up + p_lo)/2

        print(alpha_mid)

    alpha_hat = alpha_lo
    p_at_alpha_hat = p_mid


# create the grid
# ---

# we can't do a grid to infinity, so I've chosen 100 somewhat arbitrarily to produce a nice-looking graph
alpha_grid_max = alpha_hat if not np.isinf(alpha_hat) else 100

if alpha_0 < alpha_1:
    alphaV = [0, alpha_0] + list(np.linspace(alpha_0, alpha_1, ngrid//2)[1:]) + list(np.linspace(alpha_1, alpha_grid_max, ngrid)[1:])
else:
    alphaV = [0, alpha_1] + list(np.linspace(alpha_1, alpha_0, ngrid//2)[1:]) + list(np.linspace(alpha_0, alpha_grid_max, ngrid)[1:])

p_stableV = [1, 1]
p_unstableV = [0, 0]


# fill out the rest of the isoclines
# ---

print('-----')
print('finding isoclines')

for alpha in alphaV[2:]:

    # stable isocline

    if alpha <= alpha_1:

        p_stableV.append(1)

    else:

        deltap_lo = calc_delta_p_members_attract(parD, lm, partns, M, alpha, p_at_alpha_hat)
        deltap_hi = calc_delta_p_members_attract(parD, lm, partns, M, alpha, 1-1e-6)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_members_attract(parD, lm, partns, M, alpha, p), p_at_alpha_hat, 1-1e-6)
        else:
            ps = np.nan

        p_stableV.append(ps)


    # unstable isocline

    if alpha <= alpha_0:

        p_unstableV.append(0)

    else:

        deltap_lo = calc_delta_p_members_attract(parD, lm, partns, M, alpha, 1e-6)
        deltap_hi = calc_delta_p_members_attract(parD, lm, partns, M, alpha, p_at_alpha_hat)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_members_attract(parD, lm, partns, M, alpha, p), 1e-6, p_at_alpha_hat)
        else:
            ps = np.nan

        p_unstableV.append(ps)


# write results to csv
# ---

# write the alphaV and isoclines 
df_out = pd.DataFrame(list(zip(alphaV, p_stableV, p_unstableV)), columns=['alpha', 'p_stable', 'p_unstable'])
df_out.to_csv(res_dir + 'isocline' + suffix + '.csv', index=False)

if False: # NOTE
    # append the critical points to file
    dfD = [{
            'suffix': suffix,
            'n': n,
            'tau': tau,
            'W': W,
            'X': X,
            'Y': Y,
            'Z': Z,
            'alpha_0': alpha_0,
            'alpha_1': alpha_1,
            'alpha_hat': alpha_hat,
            'p_at_alpha_hat': p_at_alpha_hat,
            }]
    df_out = pd.DataFrame.from_records(dfD)
    fname = res_dir + 'isocline_critpts.csv'
    if not os.path.isfile(fname):
        # write with headers
        df_out.to_csv(fname, mode='w', header=True, index=False)
    else:
        # append
        df_out.to_csv(fname, mode='a', header=False, index=False)

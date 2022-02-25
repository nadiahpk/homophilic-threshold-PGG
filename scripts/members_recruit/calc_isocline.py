# create a detailed example

import os
from math import factorial
from scipy.special import binom
import numpy as np
from scipy.optimize import bisect
import pandas as pd

import sys
sys.path.append('../../functions/')

from my_functions import read_matrix_M, get_FV_members_recruit, calc_p_peak_random_group, calc_delta_p_members_recruit, q0_root_eqn_members_recruit, q1_root_eqn_members_recruit

# ----------------------

# user parameters
# ---

# how many increments on the x-axis for plotting
ngrid = 50

# which default parameter values to use from the file default_params.csv
suffix = '_6'
suffix = '_7'
suffix = '_5'

res_dir = '../../results/members_recruit/'

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



# find q_0, the maximum mistake probability where collectivists can still invade
# ---

# first, verify a high and low
q_lo = 1e-6
q_hi = 1-1e-6

if np.sign(q0_root_eqn_members_recruit(parD, lm, partns, M, q_lo)) != np.sign(q0_root_eqn_members_recruit(parD, lm, partns, M, q_hi)):

    q_0 = bisect(lambda q: q0_root_eqn_members_recruit(parD, lm, partns, M, q), q_lo, q_hi)

else:

    q_0 = np.nan

print('q_0 = ' + str(q_0))


# find q_1, the minimum mistake probability where individualists can invade
# ---

# first, verify a high and low
q_lo = 1e-6
q_hi = 1-1e-6

if np.sign(q1_root_eqn_members_recruit(parD, lm, partns, M, q_lo)) != np.sign(q1_root_eqn_members_recruit(parD, lm, partns, M, q_hi)):

    q_1 = bisect(lambda q: q1_root_eqn_members_recruit(parD, lm, partns, M, q), q_lo, q_hi)

else:

    q_1 = np.nan

print('q_1 = ' + str(q_1))


# find p_peak (provides bound for searching for isoclines)
# ---

p_peak = calc_p_peak_random_group(parD)


# find the maximum q where there are two isoclines
# ---

# if Delta p at p_peak < 0, then there are zero interior equilibria when q=1;
# otherwise, there are two
delta_p_at_q1 = calc_delta_p_members_recruit(parD, lm, partns, M, 1, p_peak)

if delta_p_at_q1 > 0:

    q_hat = 1
    p_at_q_hat = p_peak

else:

    # find where the transition from 2 to 0 equilibria is
    # ---

    print('-----')
    print('finding where the transition from 2 to 0 equilibria is')

    q_hi = 1
    q_lo = q_1 if q_1 > q_0 else q_0
    q_mid = (q_lo + q_hi)/2
    p_mid = p_peak

    while abs(q_lo-q_hi) > 1e-4:

        # try to find the upper and lower isoclines at q_mid

        try:
            p_up = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q_mid, p), p_mid, 1-1e-6)
        except:
            p_up = np.nan

        try:
            p_lo = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q_mid, p), 1e-6, p_mid)
        except:
            p_lo = np.nan

        if np.isnan(p_lo) or np.isnan(p_lo):
            q_hi = q_mid                # search to the left
            q_mid = (q_lo + q_mid)/2
        else:
            q_lo = q_mid                # search to the right
            q_mid = (q_hi + q_mid)/2
            p_mid = (p_up + p_lo)/2

        print(q_mid)

    q_hat = q_lo
    p_at_q_hat = p_mid


# create the grid
# ---

if q_0 < q_1:
    qV = [0, q_0] + list(np.linspace(q_0, q_1, ngrid//2)[1:]) + list(np.linspace(q_1, q_hat, ngrid)[1:])
else:
    qV = [0, q_1] + list(np.linspace(q_1, q_0, ngrid//2)[1:]) + list(np.linspace(q_0, q_hat, ngrid)[1:])

p_stableV = [1, 1]
p_unstableV = [0, 0]


# fill out the rest of the isoclines
# ---

print('-----')
print('finding isoclines')

for q in qV[2:]:

    # stable isocline

    if q <= q_1:

        p_stableV.append(1)

    else:

        deltap_lo = calc_delta_p_members_recruit(parD, lm, partns, M, q, p_at_q_hat)
        deltap_hi = calc_delta_p_members_recruit(parD, lm, partns, M, q, 1-1e-6)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q, p), p_at_q_hat, 1-1e-6)
        else:
            ps = np.nan

        p_stableV.append(ps)


    # unstable isocline

    if q <= q_0:

        p_unstableV.append(0)

    else:

        deltap_lo = calc_delta_p_members_recruit(parD, lm, partns, M, q, 1e-6)
        deltap_hi = calc_delta_p_members_recruit(parD, lm, partns, M, q, p_at_q_hat)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q, p), 1e-6, p_at_q_hat)
        else:
            ps = np.nan

        p_unstableV.append(ps)


# write results to csv
# ---

# write the qV and isoclines 
df_out = pd.DataFrame(list(zip(qV, p_stableV, p_unstableV)), columns=['q', 'p_stable', 'p_unstable'])
df_out.to_csv(res_dir + 'isocline' + suffix + '.csv', index=False)

# append the critical points to file
dfD = [{
        'suffix': suffix,
        'n': n,
        'tau': tau,
        'W': W,
        'X': X,
        'Y': Y,
        'Z': Z,
        'q_0': q_0,
        'q_1': q_1,
        'q_hat': q_hat,
        'p_at_q_hat': p_at_q_hat,
        }]
df_out = pd.DataFrame.from_records(dfD)
fname = res_dir + 'isocline_critpts.csv'
if not os.path.isfile(fname):
    # write with headers
    df_out.to_csv(fname, mode='w', header=True, index=False)
else:
    # append
    df_out.to_csv(fname, mode='a', header=False, index=False)

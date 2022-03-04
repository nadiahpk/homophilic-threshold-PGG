# create a detailed example

import os
from math import factorial
from scipy.special import binom
import numpy as np
from scipy.optimize import bisect
import pandas as pd

import sys
sys.path.append('../../functions/')

from my_functions import read_matrix_M, calc_delta_p_members_recruit
#, get_FV_members_recruit, calc_p_peak_random_group, q0_root_eqn_members_recruit, q1_root_eqn_members_recruit

# ----------------------

# threshold game payoff functions

def a_sig(parD, k):

    # unpack parameter values
    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']
    s = parD['s']

    # calculate payoff 
    l = lambda k: 1/(1+np.exp(s*(tau-1.5-k)/(n-1)))
    a_k = X + (W-X)*(l(k)-l(0))/(l(n-1)-l(0))

    return a_k

def b_sig(parD, k):

    # unpack parameter values
    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']
    s = parD['s']

    # calculate payoff 
    l = lambda k: 1/(1+np.exp(s*(tau-0.5-k)/(n-1)))
    b_k = Z + (Y-Z)*(l(k)-l(0))/(l(n-1)-l(0))

    return b_k

# ----------------------

# user parameters
# ---

s = 25
s = 5

# how many increments on the x-axis for plotting
ngrid = 50

# which default parameter values to use from the file default_params.csv
suffix = '_6'
suffix = '_7'
suffix = '_5'

res_dir = '../../results/sigmoid/'


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

# add s
parD['s'] = s


# get previously saved lm, partns, and M
# ---

lm, partns, M_num, M_den = read_matrix_M('../../results/matrix_M/matrix_M' + str(n) + '.csv')
M = M_num / M_den # the matrix for converting F to theta


# define a(k) and b(k) for this parameter set
# ---

a_fnc = lambda k: a_sig(parD, k)
b_fnc = lambda k: b_sig(parD, k)


# find q_0, the maximum mistake probability where cooperators can still invade
# ---

q_lo = 1e-6
q_hi = 1-1e-6

p_lo = 1e-6 # approximate invasion of cooperators
q_0_fnc = lambda q: calc_delta_p_members_recruit(parD, lm, partns, M, q, p_lo, a=a_fnc, b=b_fnc)

# find where q_0_func changes sign
q_0 = bisect(q_0_fnc, q_lo, q_hi)

print('q_0 = ' + str(q_0))


# find q_1, the minimum mistake probability where individualists can invade
# ---

q_lo = 1e-6
q_hi = 1-1e-6

p_hi = 1-1e-6 # approximate invasion of defectors
q_1_fnc = lambda q: calc_delta_p_members_recruit(parD, lm, partns, M, q, p_hi, a=a_fnc, b=b_fnc)

# find where q_0_func changes sign
q_1 = bisect(q_1_fnc, q_lo, q_hi)
print('q_1 = ' + str(q_1))


# find p_peak (provides bound for searching for isoclines)
# ---

# this is useful for separating the regions for searching for unstable and stable equilibria
# the peak occurs when the derivative of the gain function = 0

dV = [ a_sig(parD, k) - b_sig(parD, k) for k in range(n) ]

# g(p) = sum( binom(n-1, k) * p**k * (1-p)**(n-1-k) * dV[k] for k in range(n) )
# therefore, g'(p) = 
dgdp = lambda p: sum( binom(n-1, k) * k*p**(k-1) * (1-p)**(n-1-k) * dV[k] for k in range(n) ) \
                 - sum( binom(n-1, k) * p**k * (n-k-1)*(1-p)**(n-2-k) * dV[k] for k in range(n) )

# solve numerically for g'(p) = 0
p_lo = 1e-6
p_hi = 1-1e-6
p_peak = bisect( dgdp, p_lo, p_hi )


# find the maximum q where there are two isoclines
# ---

# if Delta p at p_peak < 0, then there are zero interior equilibria when q=1;
# otherwise, there are two
delta_p_at_q1 = calc_delta_p_members_recruit(parD, lm, partns, M, 1, p_peak, a=a_fnc, b=b_fnc)

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
            p_up = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q_mid, p, a=a_fnc, b=b_fnc), p_mid, 1-1e-6)
        except:
            p_up = np.nan

        try:
            p_lo = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q_mid, p, a=a_fnc, b=b_fnc), 1e-6, p_mid)
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

        deltap_lo = calc_delta_p_members_recruit(parD, lm, partns, M, q, p_at_q_hat, a=a_fnc, b=b_fnc)
        deltap_hi = calc_delta_p_members_recruit(parD, lm, partns, M, q, 1-1e-6, a=a_fnc, b=b_fnc)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q, p, a=a_fnc, b=b_fnc), p_at_q_hat, 1-1e-6)
        else:
            ps = np.nan

        p_stableV.append(ps)


    # unstable isocline

    if q <= q_0:

        p_unstableV.append(0)

    else:

        deltap_lo = calc_delta_p_members_recruit(parD, lm, partns, M, q, 1e-6, a=a_fnc, b=b_fnc)
        deltap_hi = calc_delta_p_members_recruit(parD, lm, partns, M, q, p_at_q_hat, a=a_fnc, b=b_fnc)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q, p, a=a_fnc, b=b_fnc), 1e-6, p_at_q_hat)
        else:
            ps = np.nan

        p_unstableV.append(ps)


# write results to csv
# ---

# write the qV and isoclines 
df_out = pd.DataFrame(list(zip(qV, p_stableV, p_unstableV)), columns=['q', 'p_stable', 'p_unstable'])
df_out.to_csv(res_dir + 'isocline' + suffix + '_s_' + str(s) + '.csv', index=False)

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
fname = res_dir + 'isocline_critpts_s_' + str(s) + '.csv'
if not os.path.isfile(fname):
    # write with headers
    df_out.to_csv(fname, mode='w', header=True, index=False)
else:
    # append
    df_out.to_csv(fname, mode='a', header=False, index=False)

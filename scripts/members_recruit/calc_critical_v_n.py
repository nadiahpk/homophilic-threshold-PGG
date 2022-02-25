# plot the critical q values against n

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

res_dir = '../../results/members_recruit/'

# how many increments on the x-axis for plotting
ngrid = 50

nV =   [4, 6, 8, 10, 12, 14, 16, 18]
tauV = [2, 3, 4,  5,  6,  7,  8,  9]

# which default parameter values to use from the file default_params.csv
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
q0V = list()
q1V = list()
qhV = list()   # store q_star values

for n, tau in zip(nV, tauV):

    print('doing n = ' + str(n) + ', tau = ' + str(tau))

    # update parameter values
    parD['n'] = n
    parD['tau'] = tau

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

    # store
    q0V.append(q_0)


    # find q_1, the minimum mistake probability where individualists can invade
    # ---

    # first, verify a high and low
    q_lo = 1e-6
    q_hi = 1-1e-6

    if np.sign(q1_root_eqn_members_recruit(parD, lm, partns, M, q_lo)) != np.sign(q1_root_eqn_members_recruit(parD, lm, partns, M, q_hi)):

        q_1 = bisect(lambda q: q1_root_eqn_members_recruit(parD, lm, partns, M, q), q_lo, q_hi)

    else:

        q_1 = np.nan

    # store
    q1V.append(q_1)


    # find q_hat, if it exists
    # ---

    # find p_peak (provides bound for searching for isoclines)
    p_peak = calc_p_peak_random_group(parD)


    # find the maximum q where there are two isoclines

    # if Delta p at p_peak < 0, then there are zero interior equilibria when q=1;
    # otherwise, there are two
    delta_p_at_q1 = calc_delta_p_members_recruit(parD, lm, partns, M, 1, p_peak)

    if delta_p_at_q1 > 0:

        q_max = 1
        p_sing = p_peak

    else:

        # find where the transition from 2 to 0 equilibria is
        # ---

        q_hi = 1
        q_lo = q_1 if q_1 > q_0 else q_0
        q_mid = (q_lo + q_hi)/2
        p_mid = p_peak

        while abs(q_lo-q_hi) > 1e-4:

            # try to find the upper and lower isoclines at q_mid

            try:
                p_up = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q_mid, p), 
                        p_mid, 1-1e-6)
            except:
                p_up = np.nan

            try:
                p_lo = bisect(lambda p: calc_delta_p_members_recruit(parD, lm, partns, M, q_mid, p), 
                        1e-6, p_mid)
            except:
                p_lo = np.nan

            if np.isnan(p_lo) or np.isnan(p_lo):
                q_hi = q_mid                # search to the left
                q_mid = (q_lo + q_mid)/2
            else:
                q_lo = q_mid                # search to the right
                q_mid = (q_hi + q_mid)/2
                p_mid = (p_up + p_lo)/2

        q_max = q_lo
        p_sing = p_mid

    # store
    qhV.append(q_max)


# write isocline data to file
# ---

# nV, q0V, q1V, qhV
df_out = pd.DataFrame(list(zip(nV, q0V, q1V, qhV)), columns=['nV', 'q0V', 'q1V', 'qhV'])
fname = res_dir + 'critical_v_n' + suffix + '.csv'
df_out.to_csv(fname, index=False)

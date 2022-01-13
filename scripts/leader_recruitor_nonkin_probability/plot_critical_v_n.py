# plot the critical q values against n

from math import factorial
from scipy.special import binom
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
import pandas as pd

import sys
sys.path.append('../../functions/')

from my_functions import read_matrix_M, calc_p_peak_random_group, calc_delta_p_leader_probability, q0_root_eqn_leader_probability, q1_root_eqn_leader_probability

# ----------------------

# user parameters
# ---

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

    if np.sign(q0_root_eqn_leader_probability(parD, lm, partns, M, q_lo)) != np.sign(q0_root_eqn_leader_probability(parD, lm, partns, M, q_hi)):

        q_0 = bisect(lambda q: q0_root_eqn_leader_probability(parD, lm, partns, M, q), q_lo, q_hi)

    else:

        q_0 = np.nan

    # store
    q0V.append(q_0)


    # find q_1, the minimum mistake probability where individualists can invade
    # ---

    # first, verify a high and low
    q_lo = 1e-6
    q_hi = 1-1e-6

    if np.sign(q1_root_eqn_leader_probability(parD, lm, partns, M, q_lo)) != np.sign(q1_root_eqn_leader_probability(parD, lm, partns, M, q_hi)):

        q_1 = bisect(lambda q: q1_root_eqn_leader_probability(parD, lm, partns, M, q), q_lo, q_hi)

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
    delta_p_at_q1 = calc_delta_p_leader_probability(parD, lm, partns, M, 1, p_peak)

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
                p_up = bisect(lambda p: calc_delta_p_leader_probability(parD, lm, partns, M, q_mid, p), 
                        p_mid, 1-1e-6)
            except:
                p_up = np.nan

            try:
                p_lo = bisect(lambda p: calc_delta_p_leader_probability(parD, lm, partns, M, q_mid, p), 
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


# plot
# ---

# where cooperators can invade
plt.plot(nV, q0V, '-o', color='blue')#, label=r'$q_0$')
plt.fill_between(nV, q0V, alpha=0.25, color='blue', label='cooperators can invade')

# where defectors can invade
plt.plot(nV, q1V, '-o', color='red')#, label=r'$q_1$')
plt.fill_between(nV, q1V, 1, alpha=0.25, color='red', label='defectors can invade')

# where cooperators cannot persist
plt.plot(nV, qhV, '-o', color='black')#, label=r'$\hat{q}$')
qvV = [ max([q0,qh]) for q0, qh in zip(q0V, qhV) ]
plt.fill_between(nV, qvV, 1, alpha=0.25, color='black', label='cooperators cannot persist')

# write as text each of the lines' identity

if q0V[-1] > q1V[-1]:
    q0va = 'bottom'
    q1va = 'top'
else:
    q0va = 'top'
    q1va = 'bottom'

plt.annotate(r'  $q_0$', (nV[-1], q0V[-1]), ha='left', va=q0va, fontsize='x-large', color='blue')
plt.annotate(r'  $q_1$', (nV[-1], q1V[-1]), ha='left', va=q1va, fontsize='x-large', color='red')
plt.annotate(r'  $\hat{q}$', (nV[-1], qhV[-1]), ha='left', va='bottom', fontsize='x-large', color='black')

plt.xlabel(r'group size $n$', fontsize='x-large')
plt.ylabel(r'probability to recruit nonkin $q$', fontsize='x-large')
plt.ylim((-0.02, 1.02))
plt.xlim((nV[0]-0.25, nV[-1]+0.25))
#plt.legend(loc='lower center', framealpha=1, fontsize='x-small', ncol=3)
plt.tight_layout()
plt.savefig('../../results/leader_recruitor_nonkin_probability/critical_v_n' + suffix + '.pdf')
plt.close()

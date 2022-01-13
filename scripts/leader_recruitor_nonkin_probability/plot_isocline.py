# create a detailed example

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

# which default parameter values to use from the file default_params.csv
suffix = '_5'
suffix = '_7'


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

if np.sign(q0_root_eqn_leader_probability(parD, lm, partns, M, q_lo)) != np.sign(q0_root_eqn_leader_probability(parD, lm, partns, M, q_hi)):

    q_0 = bisect(lambda q: q0_root_eqn_leader_probability(parD, lm, partns, M, q), q_lo, q_hi)

else:

    q_0 = np.nan

print('q_0 = ' + str(q_0))


# find q_1, the minimum mistake probability where individualists can invade
# ---

# first, verify a high and low
q_lo = 1e-6
q_hi = 1-1e-6

if np.sign(q1_root_eqn_leader_probability(parD, lm, partns, M, q_lo)) != np.sign(q1_root_eqn_leader_probability(parD, lm, partns, M, q_hi)):

    q_1 = bisect(lambda q: q1_root_eqn_leader_probability(parD, lm, partns, M, q), q_lo, q_hi)

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
delta_p_at_q1 = calc_delta_p_leader_probability(parD, lm, partns, M, 1, p_peak)

if delta_p_at_q1 > 0:

    q_max = 1
    p_sing = p_peak

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
            p_up = bisect(lambda p: calc_delta_p_leader_probability(parD, lm, partns, M, q_mid, p), p_mid, 1-1e-6)
        except:
            p_up = np.nan

        try:
            p_lo = bisect(lambda p: calc_delta_p_leader_probability(parD, lm, partns, M, q_mid, p), 1e-6, p_mid)
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

    q_max = q_lo
    p_sing = p_mid


# create the grid
# ---

if q_0 < q_1:
    qV = [0, q_0] + list(np.linspace(q_0, q_1, ngrid//2)[1:]) + list(np.linspace(q_1, q_max, ngrid)[1:])
else:
    qV = [0, q_1] + list(np.linspace(q_1, q_0, ngrid//2)[1:]) + list(np.linspace(q_0, q_max, ngrid)[1:])

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

        deltap_lo = calc_delta_p_leader_probability(parD, lm, partns, M, q, p_sing)
        deltap_hi = calc_delta_p_leader_probability(parD, lm, partns, M, q, 1-1e-6)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_leader_probability(parD, lm, partns, M, q, p), p_sing, 1-1e-6)
        else:
            ps = np.nan

        p_stableV.append(ps)


    # unstable isocline

    if q <= q_0:

        p_unstableV.append(0)

    else:

        deltap_lo = calc_delta_p_leader_probability(parD, lm, partns, M, q, 1e-6)
        deltap_hi = calc_delta_p_leader_probability(parD, lm, partns, M, q, p_sing)

        if np.sign(deltap_lo) != np.sign(deltap_hi):
            ps = bisect(lambda p: calc_delta_p_leader_probability(parD, lm, partns, M, q, p), 1e-6, p_sing)
        else:
            ps = np.nan

        p_unstableV.append(ps)


# plot
# ---

# plot the isoclines

plt.plot(qV, p_stableV,   color='black', lw=3, label='stable')
plt.plot(qV, p_unstableV, color='black', ls='dashed', label='unstable')

# plot the p=0 and p=1 continuation of the isoclines

plt.plot([q_0, 1], [0, 0],  color='black', lw=3)
plt.plot([q_1, 1], [1, 1],  color='black', ls='dashed')

# plot the point where the transition from 2 to 0 interior equilibria occurs (if needed)
if delta_p_at_q1 < 0:
    plt.scatter([q_max], [p_sing], color='black')
    plt.plot([q_max, q_max], [-0.1, 1], color='black', alpha=0.5, lw=1)
    plt.text(q_max, 0, r'$\hat{q}$ ' + '\n', ha='right', va='bottom', fontsize='large')

# plot where the p isocline intersects x axis

plt.plot([q_0, q_0], [-0.1, 1], color='black', alpha=0.5, lw=1)
plt.text(q_0, 0, r'$q_0$' + '\n', ha='right', va='bottom', fontsize='large')

plt.plot([q_1, q_1], [-0.1, 1], color='black', alpha=0.5, lw=1)
plt.text(q_1, 1, '\n' + r'$q_1$', ha='right', va='top', fontsize='large')

# colour regions

plt_ymin = -0.05
plt_ymax = 1.12
plt_yspan = plt_ymax - plt_ymin

opacity = 0.25
plt.axvspan(0, q_0, alpha=opacity, color='blue',
        ymin=(0-plt_ymin)/plt_yspan, ymax=(1-plt_ymin)/plt_yspan)
plt.axvspan(q_1, 1, alpha=opacity, color='red',
        ymin=(0-plt_ymin)/plt_yspan, ymax=(1-plt_ymin)/plt_yspan)
plt.axvspan(q_max, 1, alpha=opacity, color='black',
        ymin=(0-plt_ymin)/plt_yspan, ymax=(1-plt_ymin)/plt_yspan)


# decorate the plot

# explain the x-axis
plt.text(1, -0.12, r'$\uparrow$', fontsize='small', ha='center', va='top')
plt.text(1.01, -0.12, '\nrandom groups', fontsize='small', ha='right', va='top')
plt.text(0, -0.12, r'$\uparrow$', fontsize='small', ha='center', va='top')
plt.text(-0.01, -0.12, '\nperfect homophily', fontsize='small', ha='left', va='top')

plt.xlabel(r'probability to recruit nonkin $q$', fontsize='x-large')
plt.ylabel(r'proportion of Cooperators $p$', fontsize='x-large')
plt.xlim((0, 1))
plt.ylim((plt_ymin, plt_ymax))
plt.legend(loc='upper center', framealpha=1, ncol=2)

plt.tight_layout()
plt.savefig('../../results/leader_recruitor_nonkin_probability/isocline' + suffix + '.pdf')
plt.close()

# plot multiple isoclines for different group sizes n

from math import factorial
from scipy.special import binom
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect
import pandas as pd

import sys
sys.path.append('../../functions/')

from my_functions import read_matrix_M, get_FV_random_probability, calc_p_peak_random_group, calc_delta_p_random_probability, q0_root_eqn_random_probability, q1_root_eqn_random_probability

# ----------------------

# user parameters
# ---

ngrid = 50 # how many increments on the x-axis for plotting

# parameter values common to all the isoclines we'll plot
common_pars = {
        'n': 8,     # group size held constant
        'W':  2,    # cooperator payoff if threshold met
        'X': -1,    # cooperator payoff if threshold not met
        'Y':  3,    # defector payoff if threshold met
        'Z':  0,    # defector payoff if threshold not met
        }

# parameter values specific to each isocline
specific_pars = [
        {'tau': 7, 'colour': 'brown'},
        {'tau': 6, 'colour': 'magenta'},
        {'tau': 5, 'colour': 'red'},
        {'tau': 4, 'colour': 'blue'},
        {'tau': 3, 'colour': 'orange'},
        {'tau': 2, 'colour': 'black'}
        ]


# calculate each isocline corresponding to specific_pars
# ---

# each isocline is calculated in a similar way to plot_isocline.py

# a place to store info for each plot
qM = list()
p_stableM = list()
p_unstableM = list()
q_0V = list()
q_1V = list()
q_maxV = list()
p_singV = list()

for specific_par in specific_pars:

    # create parD, which is used for various functions
    parD = { par_name: par_val for par_name, par_val in common_pars.items() } # n, W, X, Y, Z
    parD['tau'] = specific_par['tau']

    # unpack needed parameter values
    n = parD['n']
    tau = parD['tau']

    print('-----')
    print('doing tau = ' + str(tau))


    # get previously saved lm, partns, and M
    # ---

    lm, partns, M_num, M_den = read_matrix_M('../../results/matrix_M/matrix_M' + str(n) + '.csv')
    M = M_num / M_den # the matrix for converting F to theta


    # find q_0, the maximum mistake probability where collectivists can still invade
    # ---

    # first, verify a high and low
    q_lo = 1e-6
    q_hi = 1-1e-6

    if np.sign(q0_root_eqn_random_probability(parD, lm, partns, M, q_lo)) != np.sign(q0_root_eqn_random_probability(parD, lm, partns, M, q_hi)):

        q_0 = bisect(lambda q: q0_root_eqn_random_probability(parD, lm, partns, M, q), q_lo, q_hi)

    else:

        q_0 = np.nan

    # store
    q_0V.append(q_0)


    # find q_1, the minimum mistake probability where individualists can invade
    # ---

    # first, verify a high and low
    q_lo = 1e-6
    q_hi = 1-1e-6

    if np.sign(q1_root_eqn_random_probability(parD, lm, partns, M, q_lo)) != np.sign(q1_root_eqn_random_probability(parD, lm, partns, M, q_hi)):

        q_1 = bisect(lambda q: q1_root_eqn_random_probability(parD, lm, partns, M, q), q_lo, q_hi)

    else:

        q_1 = np.nan

    # store
    q_1V.append(q_1)


    # find p_peak (provides bound for searching for isoclines)
    # ---

    p_peak = calc_p_peak_random_group(parD)


    # find the maximum q where there are two isoclines
    # ---

    # if Delta p at p_peak < 0, then there are zero interior equilibria when q=1;
    # otherwise, there are two
    delta_p_at_q1 = calc_delta_p_random_probability(parD, lm, partns, M, 1, p_peak)

    if delta_p_at_q1 > 0:

        q_max = 1
        p_sing = p_peak

    else:

        # find where the transition from 2 to 0 equilibria is
        # ---

        print('finding where the transition from 2 to 0 equilibria is')

        q_hi = 1
        q_lo = q_1 if q_1 > q_0 else q_0
        q_mid = (q_lo + q_hi)/2
        p_mid = p_peak

        while abs(q_lo-q_hi) > 1e-4:

            # try to find the upper and lower isoclines at q_mid

            try:
                p_up = bisect(lambda p: calc_delta_p_random_probability(parD, lm, partns, M, q_mid, p), p_mid, 1-1e-6)
            except:
                p_up = np.nan

            try:
                p_lo = bisect(lambda p: calc_delta_p_random_probability(parD, lm, partns, M, q_mid, p), 1e-6, p_mid)
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

    # store
    q_maxV.append(q_max)
    p_singV.append(p_sing)


    # create the grid
    # ---

    if q_0 < q_1:
        qV = [0, q_0] + list(np.linspace(q_0, q_1, ngrid//2)[1:]) + list(np.linspace(q_1, q_max, ngrid)[1:])
    else:
        qV = [0, q_1] + list(np.linspace(q_1, q_0, ngrid//2)[1:]) + list(np.linspace(q_0, q_max, ngrid)[1:])

    qM.append(qV)


    # fill out the rest of the isoclines
    # ---

    p_stableV = [1, 1]
    p_unstableV = [0, 0]

    print('finding isoclines')

    for q in qV[2:]:

        # stable isocline

        if q <= q_1:

            p_stableV.append(1)

        else:

            deltap_lo = calc_delta_p_random_probability(parD, lm, partns, M, q, p_sing)
            deltap_hi = calc_delta_p_random_probability(parD, lm, partns, M, q, 1-1e-6)

            if np.sign(deltap_lo) != np.sign(deltap_hi):
                ps = bisect(lambda p: calc_delta_p_random_probability(parD, lm, partns, M, q, p), p_sing, 1-1e-6)
            else:
                ps = np.nan

            p_stableV.append(ps)


        # unstable isocline

        if q <= q_0:

            p_unstableV.append(0)

        else:

            deltap_lo = calc_delta_p_random_probability(parD, lm, partns, M, q, 1e-6)
            deltap_hi = calc_delta_p_random_probability(parD, lm, partns, M, q, p_sing)

            if np.sign(deltap_lo) != np.sign(deltap_hi):
                ps = bisect(lambda p: calc_delta_p_random_probability(parD, lm, partns, M, q, p), 1e-6, p_sing)
            else:
                ps = np.nan

            p_unstableV.append(ps)

    p_stableM.append(p_stableV)
    p_unstableM.append(p_unstableV)


# plot each of the isoclines
# ---

plot_lines = list()
for qV, p_stableV, p_unstableV, q_0, q_1, q_max, p_sing, s_par in zip(qM, p_stableM, p_unstableM, q_0V, q_1V, q_maxV, p_singV, specific_pars):

    tau = s_par['tau']
    colour = s_par['colour']

    # plot the isoclines

    l1, = plt.plot(qV, p_stableV,   color=colour, alpha=0.7, lw=3)
    l2, = plt.plot(qV, p_unstableV, color=colour, alpha=0.7, ls='dashed')

    # append for the legend later
    plot_lines.append([l1, l2])

    # plot the p=0 and p=1 continuation of the isoclines

    plt.plot([q_0, 1], [0, 0],  color=colour, alpha=0.7, lw=3)
    plt.plot([q_1, 1], [1, 1],  color=colour, alpha=0.7, ls='dashed')

    # plot the point where the transition from 2 to 0 interior equilibria occurs (if needed)
    if q_max < 1:
        plt.scatter([q_max], [p_sing], color=colour)


# create legend
# ---

legend1 = plt.legend(plot_lines[5], ['stable', 'unstable'], loc='upper center', framealpha=1, ncol=2)
legend2 = plt.legend([l[0] for l in plot_lines], [r'$\tau=' + str(spar['tau']) + r'$' for spar in specific_pars], loc='center left')
plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)

# decorate and save the plot
# ---

# explain the x-axis
#plt.text(1, -0.12, r'$\uparrow$', fontsize='small', ha='center', va='top')
#plt.text(1.01, -0.12, '\nrandom groups', fontsize='small', ha='right', va='top')
#plt.text(0, -0.12, r'$\uparrow$', fontsize='small', ha='center', va='top')
#plt.text(-0.01, -0.12, '\nperfect homophily', fontsize='small', ha='left', va='top')

plt.xlabel(r'probability to recruit nonkin $q$', fontsize='x-large')
plt.ylabel(r'proportion of Cooperators $p$', fontsize='x-large')
plt.xlim((0, 1))
plt.ylim((-0.05, 1.12))

plt.tight_layout()
plt.savefig('../../results/random_recruitor_nonkin_probability/isoclines_v_tau.pdf')
plt.close()

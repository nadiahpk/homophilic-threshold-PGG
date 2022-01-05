# Plot the isoclines in p vs W, X, and Y space
# example with default parameter set "_5"

import matplotlib.pyplot as plt
from scipy.special import binom
from scipy.optimize import bisect
import pandas as pd
import numpy as np


# ------------------------------------------------------------

def p_root_eqn(p, parD):
    '''
    Find the isocline by solving this equation = 0
    '''

    n = parD['n']
    tau = parD['tau']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']

    res = X - Z + (W-X) * binom(n-1,tau-1) * p**(tau-1) + sum( binom(n-1,l) * p**l * (-1)**(l-tau) * ( binom(l-1, tau-2) * (X-W) + binom(l-1, tau-1) * (Z-Y)) for l in range(tau, n) ) 

    return res

# ------------------------------------------------------------


# parameters
# ---

suffix = '_5'
par_names = ['n', 'tau', 'W', 'X', 'Y', 'Z']    # all parameter names corresponding to columns in default_params.csv
ngrid = 100


# dictionary for plot range { vary_par_name: { dictionary of plot info } }
plotD = {
        'W': { 
            'vpar_lo': None,        # lower limit tbd by the critical W value
            'vpar_hi': 4,           # upper limit
            'vpar_in': 'Y',         # where the interior equilibria appear
            'xlabel': r'Cooperator payoff if threshold met, $W$'
            },
        'X': { 
            'vpar_lo': None,
            'vpar_hi': 0,
            'vpar_in': 'Z',
            'xlabel': r'Cooperator payoff if threshold not met, $X$'
            },
        'Y': { 
            'vpar_lo': 0,
            'vpar_hi': None,
            'vpar_in': 'W',
            'xlabel': r'Defector payoff if threshold met, $Y$'
            },
        }



for vary_par_name in ['W', 'X', 'Y']:

    # read in the default and critical values
    # ---

    fname = '../../results/random_grp/find_crit' + suffix + '.csv'
    df = pd.read_csv(fname)
    df.set_index('vary_par_name', inplace=True)

    # default values
    parD_orig = { par_name: df.loc[vary_par_name][par_name] for par_name in par_names }

    # critical values
    par_crit = df.loc[vary_par_name]['vary_par_critical_value']
    p_crit = df.loc[vary_par_name]['p_critical_value']


    # the stable isocline
    # ---

    if vary_par_name == 'W' or vary_par_name == 'X':


        # go from the lower bound to vpar_int
        # ---

        vpar_lo = par_crit
        vpar_hi = parD_orig[plotD[vary_par_name]['vpar_in']]
        vary_par_vals = list(np.linspace(vpar_lo, vpar_hi, ngrid)[1:-1])

        pV = list()
        for vary_par_val in vary_par_vals:

            # update the dictionary
            parD = { par_name: par_val for par_name, par_val in parD_orig.items() }
            parD[vary_par_name] = vary_par_val

            # find the isocline point
            p = bisect( lambda p: p_root_eqn(p, parD), p_crit, 1 )
            pV.append(p)

        if vary_par_name == 'W':

            # append on the stable upper isocline
            # ---

            vary_par_vals.append(vpar_hi)
            pV.append(1)
            vary_par_vals.append(plotD[vary_par_name]['vpar_hi'])
            pV.append(1)


        # preappend the critical values
        # ---

        stable_par_vals = [par_crit] + vary_par_vals
        stable_ps = [p_crit] + pV

    elif vary_par_name == 'Y':

        # stable upper isocline from lower bound to vpar_int
        # ---

        stable_par_vals = [ plotD[vary_par_name]['vpar_lo'], parD_orig[plotD[vary_par_name]['vpar_in']] ]
        stable_ps = [1, 1]


        # go from vpar_in to upper bound
        # ---

        vpar_lo = parD_orig[plotD[vary_par_name]['vpar_in']]
        vpar_hi = par_crit
        vary_par_vals = list(np.linspace(vpar_lo, vpar_hi, ngrid)[1:-1])

        pV = list()
        for vary_par_val in vary_par_vals:

            # update the dictionary
            parD = { par_name: par_val for par_name, par_val in parD_orig.items() }
            parD[vary_par_name] = vary_par_val

            # find the isocline point
            p = bisect( lambda p: p_root_eqn(p, parD), p_crit, 1 )
            pV.append(p)

        stable_par_vals = stable_par_vals + vary_par_vals
        stable_ps = stable_ps + pV

        # append the critical values
        # ---

        stable_par_vals.append(par_crit)
        stable_ps.append(p_crit)


    # the unstable isocline
    # ---

    if vary_par_name == 'W' or vary_par_name == 'X':

        vpar_lo = par_crit
        vpar_hi = plotD[vary_par_name]['vpar_hi']
        vary_par_vals = list(np.linspace(vpar_lo, vpar_hi, ngrid)[1:])

        pV = list()
        for vary_par_val in vary_par_vals:

            # update the dictionary
            parD = { par_name: par_val for par_name, par_val in parD_orig.items() }
            parD[vary_par_name] = vary_par_val

            # find the isocline point
            p = bisect( lambda p: p_root_eqn(p, parD), 0, p_crit )
            pV.append(p)

        # preappend critical
        unstable_par_vals = [par_crit] + vary_par_vals
        unstable_ps = [p_crit] + pV

    else:

        vpar_lo = plotD[vary_par_name]['vpar_lo']
        vpar_hi = par_crit
        vary_par_vals = list(np.linspace(vpar_lo, vpar_hi, ngrid)[:-1])

        pV = list()
        for vary_par_val in vary_par_vals:

            # update the dictionary
            parD = { par_name: par_val for par_name, par_val in parD_orig.items() }
            parD[vary_par_name] = vary_par_val

            # find the isocline point
            p = bisect( lambda p: p_root_eqn(p, parD), 0, p_crit )
            pV.append(p)

        # postappend critical
        unstable_par_vals = vary_par_vals + [par_crit] 
        unstable_ps = pV + [p_crit] 


    # plot isoclines
    # ---

    fname = '../../results/random_grp/isoclines_' + vary_par_name + suffix + '.pdf'

    # get x-range of plot
    if plotD[vary_par_name]['vpar_lo'] is None:
        vpar_lo = np.floor(par_crit)
    else:
        vpar_lo = plotD[vary_par_name]['vpar_lo']

    if plotD[vary_par_name]['vpar_hi'] is None:
        vpar_hi = np.ceil(par_crit)
    else:
        vpar_hi = plotD[vary_par_name]['vpar_hi']

    # unstable isocline at p=1
    plt.plot([vpar_lo, vpar_hi], [1, 1], color='black', ls='dashed')

    # stable isocline at p=0
    plt.plot([vpar_lo, vpar_hi], [0, 0], color='black', lw=3)

    # critical point
    plt.scatter([par_crit], [p_crit], color='black')

    # where two interior equilibria appear
    xloc = parD_orig[plotD[vary_par_name]['vpar_in']]
    pretty = r' $' + plotD[vary_par_name]['vpar_in'] + r'$'
    plt.plot([xloc, xloc], [-0.1, 1], color='black', alpha=0.5, lw=1)
    plt.text(xloc, 0.05, pretty, ha='left', va='bottom')

    # stable and unstable isocline
    plt.plot(stable_par_vals, stable_ps, color='black', lw=3, label='stable')
    plt.plot(unstable_par_vals, unstable_ps, color='black', ls='dashed', label='unstable')

    plt.xlabel(plotD[vary_par_name]['xlabel'], fontsize='x-large')
    plt.ylabel(r'proportion of Cooperators, $p$', fontsize='x-large')

    if vary_par_name == 'X':
        plt.xlim((vpar_lo, vpar_hi+0.05))
    else:
        plt.xlim((vpar_lo, vpar_hi))

    plt.ylim((-0.05, 1.12))


    # Hisashi wants some arrows
    # ---

    def make_arrow(x, y, dy):
        plt.annotate('', xy=(x, y), xytext=(x, y+dy), arrowprops=dict(arrowstyle='->'))

    if vary_par_name == 'W':

        make_arrow(2.50, 0.10, +0.1)
        make_arrow(2.50, 0.62, -0.1)
        make_arrow(2.50, 0.78, +0.1)

        make_arrow(3.50, 0.10, +0.1)
        make_arrow(3.50, 0.80, -0.1)

        make_arrow(2.07, 0.54, +0.1)

    if vary_par_name == 'Y':

        make_arrow(2.25, 0.10, +0.1)
        make_arrow(2.25, 0.67, -0.1)
        make_arrow(2.25, 0.84, +0.1)

        make_arrow(1.00, 0.10, +0.1)
        make_arrow(1.00, 0.80, -0.1)

        make_arrow(2.82, 0.55, +0.1)

    if vary_par_name == 'X':

        make_arrow(-0.3, 0.75, +0.1)
        make_arrow(-0.3, 0.55, -0.1)
        make_arrow(-0.3, 0.10, +0.1)

        make_arrow(-0.9, 0.5, +0.1)

    plt.legend(loc='upper center', framealpha=1, ncol=2)
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()

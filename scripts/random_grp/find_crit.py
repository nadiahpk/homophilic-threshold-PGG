# The critical values of W, X, Y, or Z are the values at which the 
# evolutionary dynamics transitions from having 0 interior equilibria 
# to having 2 equilibria. We take one parameter (e.g., W) while holding
# the other parameters (e.g., X, Y, Z) constant. We first find where 
# p_peak is, which is the maximum in Delta p, and then we solve for 
# when Delta p = 0.

from scipy.special import binom
from scipy.optimize import newton
import pandas as pd

# ------------------------------------------------------------

def calc_p_peak(tau, n, W, X, Y, Z):
    '''
    Find the value of p that maximises Delta p

    Inputs
    ---

    tau, int
        Threshold number of cooperators for public good to be produced

    n, int
        Group size

    W, float
        Cooperator payoff if threshold met

    X, float
        Cooperator payoff if threshold not met

    Y, float
        Defector payoff if threshold met

    Z, float
        Defector payoff if threshold not met


    Outputs
    ---

    p_peak, float [0,1]
        Value of p that maximises Delta p
    '''

    p_peak = (tau-1)*(W-X)/((tau-1)*(W-X)+(Y-Z)*(n-tau))

    return p_peak

# find the root of this equation to find the critical value of parameter par
def make_0(par_val, par_name, parD_orig):
    '''

    Inputs:
    ---

    par_val, float
        The value of the parameter

    par_name, string
        The name of the parameter (W, X, Y, or Z)

    parD_orig, dict {string: float}
        A dictionary of default parameter values

    Outputs:
    ---

    res, float
        The value of Delta p at the default parameter values when 
        parameter par_name has value par_val
    '''

    # make a copy and update the parameter value dictionary
    parD = { par_name: par_val for par_name, par_val in parD_orig.items() }
    parD[par_name] = par_val

    # unpack parameter values
    tau = parD['tau']
    n = parD['n']
    W = parD['W']
    X = parD['X']
    Y = parD['Y']
    Z = parD['Z']

    # find p_peak: (tau-1)(W-X)/[(tau-1)(W-X)+(Y-Z)(n-tau)]
    p_peak = calc_p_peak(tau, n, W, X, Y, Z)

    # Delta p given p_peak equals 0 at the critical par value
    res = X - Z + (W-X) * binom(n-1,tau-1) * p_peak**(tau-1) + sum( binom(n-1,l) * p_peak**l * (-1)**(l-tau) * ( binom(l-1, tau-2) * (X-W) + binom(l-1, tau-1) * (Z-Y)) for l in range(tau, n) )

    return res
# ------------------------------------------------------------

# user parameters
# ---

#suffix = '_5'                                   # which default parameter values to use from the file default_params.csv
suffix = '_6'                                   # which default parameter values to use from the file default_params.csv


# fixed parameters
# ---

par_names = ['n', 'tau', 'W', 'X', 'Y', 'Z']    # all parameter names corresponding to columns in default_params.csv
vary_par_names = ['W', 'X', 'Y']                # which parameters are we finding the critical value for


# read in the default parameter values
# ---

df = pd.read_csv('../default_params.csv')
df.set_index('suffix', inplace=True)
parD_orig = { par_name: df[par_name].loc[suffix] for par_name in par_names }


# find the critical value for each vary_par_name while holding other parameters at default values
# ---

vary_par_val_crits = list()
p_crits = list()
for vary_par_name in vary_par_names:

    # to find critical value of vary_par_name, find root of make_0()
    vary_par_val0 = parD_orig[vary_par_name] # start value is the default value
    vary_par_val_crit = newton(lambda vary_par_val: make_0(vary_par_val, vary_par_name, parD_orig), vary_par_val0)
    vary_par_val_crits.append(vary_par_val_crit)

    # critical value of p
    parD = { par_name: par_val for par_name, par_val in parD_orig.items() }
    parD[vary_par_name] = vary_par_val_crit
    tau = parD['tau']; n = parD['n']; W = parD['W']; X = parD['X']; Y = parD['Y']; Z = parD['Z'] # unpack parameter values
    p_crit = calc_p_peak(tau, n, W, X, Y, Z)
    p_crits.append(p_crit)


# write results to a file
# ---

num_rows = len(vary_par_names)
dfD = {
        'suffix': [suffix]*num_rows,
        'n': [parD_orig['n']]*num_rows,
        'tau': [parD_orig['tau']]*num_rows,
        'W': [parD_orig['W']]*num_rows,
        'X': [parD_orig['X']]*num_rows,
        'Y': [parD_orig['Y']]*num_rows,
        'Z': [parD_orig['Z']]*num_rows,
        'vary_par_name': vary_par_names,
        'vary_par_critical_value': vary_par_val_crits,
        'p_critical_value': p_crits,
        }
df = pd.DataFrame(data = dfD) # suffix, default_n, default_tau, default_W, default_X, default_Y, default_Z, par_name, par_critical_value, p_critical_value
df.to_csv('../../results/random_grp/find_crit' + suffix + '.csv', index=False)

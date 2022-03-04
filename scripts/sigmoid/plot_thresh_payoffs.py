# I need to first show readers what the payoffs look like for the two defaults

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import binom

# user parameters
# ---

# the two examples used in the main text
suffix = '_5'
suffix = '_7'

res_dir = '../../results/sigmoid/'


# threshold game payoff functions
# ---

b = lambda k: Z if k < tau else Y
a = lambda k: X if k < tau-1 else W


# read in the default parameter values
# ---

# all parameter names corresponding to columns in default_params.csv
par_names = ['n', 'tau', 'W', 'X', 'Y', 'Z']

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


# plot a_k and b_k
# ---

height = 3.8 # default 4.8
width = 1.3*height
plt.figure(figsize=(width, height)) # default

kV = list(range(n))                         # 0, ..., n-1
aV = [ a(k) for k in kV ]                   # cooperator payoff
bV = [ b(k) for k in kV ]                   # defector payoff
dV = [ ai-bi for ai, bi in zip(aV, bV) ]    # gain sequence

plt.plot(kV, aV, '-o', color='blue', label='Cooperator')
plt.plot(kV, bV, '-o', color='red',  label='Defector')
plt.xlabel(r'number of Cooperators $k$ of $n-1$', fontsize='x-large')
plt.ylabel(r'payoff', fontsize='x-large')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(res_dir + 'thresh_payoffs' + suffix + '.pdf')
plt.close()


# plot gains
# ---

height = 3.8 # default 4.8
width = 1.3*height
plt.figure(figsize=(width, height)) # default

# plot d_k

plt.plot(kV, dV, '-o', color='blue', label=r'gain sequence $d_k$')

# plot g(p)

pV = np.linspace(0, 1, 50)
gV = [ sum( binom(n-1, k) * p**k * (1-p)**(n-1-k) * dV[k] for k in range(n) ) for p in pV ]
plt.plot(pV*(n-1), gV, color='black', label=r'gain function $g(p)$')

# decorate plot

plt.xlabel(r'number of Cooperators $k$ of $n-1$', fontsize='x-large')
plt.ylabel(r'gain function and sequence', fontsize='x-large')
plt.axhline(0, color='black', alpha=0.7, lw=1)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(res_dir + 'thresh_payoffs_switch' + suffix + '.pdf')
plt.close()

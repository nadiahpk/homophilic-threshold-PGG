# put details in suffix = 5

import pandas as pd
import matplotlib.pyplot as plt

# parameters
# ---

#suffix = '_6'
suffix = '_7'
suffix = '_5'
#suffix = '_4242'
res_dir = '../../results/members_attract/'


# read in the data
# ---

fname = res_dir + 'critical_v_n' + suffix + '.csv'
df = pd.read_csv(fname)
nV = df['nV'].values
alpha0V = df['alpha0V'].values
alpha1V = df['alpha1V'].values
alphahV = df['alphahV'].values

# invert to match the new "homophily" as the variable
fnc = lambda x: 1.5**(-x)
h0V = fnc(alpha0V)
h1V = fnc(alpha1V)
hhV = fnc(alphahV)

# plot size
height = 3.7 # default 4.8
width = 1.1*height
plt.figure(figsize=(width, height)) # default

# where cooperators can invade
plt.plot(nV, h0V, '-o', color='blue')#, label=r'$h_0$')
plt.fill_between(nV, h0V, 1, alpha=0.25, color='blue', label='cooperators can invade')

# where defectors can invade
plt.plot(nV, h1V, '-o', color='red')#, label=r'$h_1$')
plt.fill_between(nV, h1V, 0, alpha=0.25, color='red', label='defectors can invade')

# where cooperators cannot persist
plt.plot(nV, hhV, '-o', color='black')#, label=r'$\hat{h}$')
hvV = [ min([h0,hh]) for h0, hh in zip(h0V, hhV) ]
plt.fill_between(nV, hvV, 0, alpha=0.25, color='black', label='cooperators cannot persist')

# write as text each of the lines' identity

if h0V[-1] > h1V[-1]:
    h0va = 'bottom'
    h1va = 'top'
else:
    h0va = 'top'
    h1va = 'bottom'

plt.annotate(r'  $h_1$ ', (nV[1], h1V[1]), ha='right', va=h1va, fontsize='x-large', color='red')
plt.annotate(r'  $h_0$', (nV[1], h0V[1]), ha='left', va=h0va, fontsize='x-large', color='blue')
plt.annotate(r'  $\hat{h}$ ', (nV[2], hhV[2]), ha='right', va='bottom', fontsize='x-large', color='black')

# finish off

plt.xlabel(r'group size $n$', fontsize='x-large')
plt.ylabel(r'homophily $h$', fontsize='x-large')
plt.ylim((-0.02, 1.02))
plt.xlim((nV[0]-0.25, nV[-1]+0.25))
plt.xticks([4, 18])
#plt.legend(loc='lower center', framealpha=1, fontsize='x-small', ncol=3)
plt.tight_layout()
plt.savefig(res_dir + 'critical_v_n' + suffix + '.pdf')
plt.close()

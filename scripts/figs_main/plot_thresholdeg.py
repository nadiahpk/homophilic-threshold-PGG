# draw an example of the threshold game for main

import matplotlib.pyplot as plt


# parameters
# ---

out_dir = '../../results/figs_main/'

# use the same parameter values as main, but hide the numbers later
n = 8
tau = 5
W = 2
X = -1
Y = 3
Z = 0


# total payoff
# ---

cV = list(range(n+1)) # total number of cooperators
totV = [ 0 if c < tau else Y*n for c in cV ]


# coop and defect payoffs
# ---

copV = [ X if c < tau else W for c in cV ]
defV = [ Z if c < tau else Y for c in cV ]


# plotting
# ---

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(3,3))


# plot total payoff

axs[0].set_ylabel('public good', labelpad=10)
axs[0].plot( cV, totV, lw=3, color='black' )
#axs[0].axhline(0, lw=1, color='black')
axs[0].tick_params(right=True, labelleft=False, labelright=True)
axs[0].set_yticks([Z, Y*n])
axs[0].set_yticklabels(['0', ''])


# plot cooperators and defectors

axs[1].set_ylabel('payoffs', labelpad=10)
axs[1].plot( cV, copV, lw=3, color='blue' )
axs[1].plot( cV, defV, lw=3, color='red' )
#axs[1].axhline(0, lw=1, color='black')
axs[1].tick_params(right=True, labelright=True, labelleft=False)
axs[1].set_yticks([W, X, Y, Z])
axs[1].set_yticklabels([r'$W$', r'$X$', r'$Y$', r'$Z$'])


# decorate

plt.xticks(cV)
#plt.xlim((-0.5, n+1.5))
plt.xlabel('number of cooperators', fontsize='large')
plt.tight_layout()
plt.savefig(out_dir + 'thresholdeg.pdf')
plt.close()

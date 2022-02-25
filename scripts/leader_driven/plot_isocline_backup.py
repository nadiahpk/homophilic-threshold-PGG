import matplotlib.pyplot as plt
import pandas as pd


# user parameters
# ---

suffix = '_7'
res_dir = '../../results/leader_driven/'


# read in the needed data
# ---

# get the critical points
fname = res_dir + 'isocline_critpts' + suffix + '.csv'
df = pd.read_csv(fname)
q_0 = df[df['suffix']==suffix].iloc[0]['q_0']
q_1 = df[df['suffix']==suffix].iloc[0]['q_1']
q_hat = df[df['suffix']==suffix].iloc[0]['q_hat']
p_at_q_hat = df[df['suffix']==suffix].iloc[0]['p_at_q_hat']

# get the isoclines
fname = res_dir + 'isocline' + suffix + '.csv'
df = pd.read_csv(fname)
qV = df['q']
p_stableV = df['p_stable']
p_unstableV = df['p_unstable']


# plot
# ---

# plot the isoclines

plt.plot(qV, p_stableV,   color='black', lw=3, label='stable')
plt.plot(qV, p_unstableV, color='black', ls='dashed', label='unstable')

# plot the p=0 and p=1 continuation of the isoclines

plt.plot([q_0, 1], [0, 0],  color='black', lw=3)
plt.plot([q_1, 1], [1, 1],  color='black', ls='dashed')

# plot the point where the transition from 2 to 0 interior equilibria occurs (if needed)
if q_hat != 1:
    plt.scatter([q_hat], [p_at_q_hat], color='black')
    plt.plot([q_hat, q_hat], [-0.1, 1], color='black', alpha=0.5, lw=1)
    plt.text(q_hat, 0, r'$\hat{q}$ ' + '\n', ha='right', va='bottom', fontsize='large')

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
plt.axvspan(q_hat, 1, alpha=opacity, color='black',
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
plt.savefig(res_dir + 'isocline' + suffix + '.pdf')
plt.close()


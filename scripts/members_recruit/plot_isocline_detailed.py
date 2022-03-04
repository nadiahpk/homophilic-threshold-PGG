# this one includes arrows to show how the isoclines attract or repel

import matplotlib.pyplot as plt
import pandas as pd


# user parameters
# ---

suffix = '_5'
suffix = '_7'
res_dir = '../../results/members_recruit/'


# read in the needed data
# ---

# get the critical points
fname = res_dir + 'isocline_critpts.csv'
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


# transform everything from q to h = 1-q
# ---

hV = 1-qV
h_0 = 1-q_0
h_1 = 1-q_1
h_hat = 1-q_hat


# plot
# ---

height = 3.8 # default 4.8
width = 1.3*height
plt.figure(figsize=(width, height)) # default

# plot the isoclines against 1-qV
qV = 1-qV

plt.plot(hV, p_stableV,   color='black', lw=3, label='stable')
plt.plot(hV, p_unstableV, color='black', ls='dashed', label='unstable')

# plot the p=0 and p=1 continuation of the isoclines

plt.plot([0, h_0], [0, 0],  color='black', lw=3)
plt.plot([0, h_1], [1, 1],  color='black', ls='dashed')

# plot the point where the transition from 2 to 0 interior equilibria occurs (if needed)
if h_hat != 0:
    plt.scatter([h_hat], [p_at_q_hat], color='black')
    plt.plot([h_hat, h_hat], [-0.1, 1], color='black', alpha=0.5, lw=1)
    #plt.text(h_hat, 0, r' $\hat{h}$ ' + '\n', ha='left', va='bottom', fontsize='x-large')

# plot where the p isocline intersects x axis

plt.plot([h_0, h_0], [-0.1, 1], color='black', alpha=0.5, lw=1)
#plt.text(h_0+0.01, 0, r'$h_0$' + '\n', ha='left', va='bottom', fontsize='x-large')

plt.plot([h_1, h_1], [-0.1, 1], color='black', alpha=0.5, lw=1)
#plt.text(h_1+0.01, 1, '\n' + r'$h_1$', ha='left', va='top', fontsize='x-large')

# colour regions

plt_ymin = -0.05
plt_ymax = 1.15
plt_yspan = plt_ymax - plt_ymin

opacity = 0.25
plt.axvspan(1, 1-q_0, alpha=opacity, color='blue',
        ymin=(0-plt_ymin)/plt_yspan, ymax=(1-plt_ymin)/plt_yspan)
plt.axvspan(1-q_1, 0, alpha=opacity, color='red',
        ymin=(0-plt_ymin)/plt_yspan, ymax=(1-plt_ymin)/plt_yspan)
plt.axvspan(1-q_hat, 0, alpha=opacity, color='black',
        ymin=(0-plt_ymin)/plt_yspan, ymax=(1-plt_ymin)/plt_yspan)


# label p_u* and p_s*

h_s = h_hat + (h_1-h_hat)/4                         # hoped position for p_s
idx_s = list(abs(hV-h_s)).index(min(abs(hV-h_s)))   # find a close value in hV
p_s = p_stableV[idx_s]
plt.text(h_s, p_s, '$p_s^*$ ', ha='right', va='bottom', fontsize='x-large')

h_u = h_hat + (h_0-h_hat)/4                         # hoped position for p_u
idx_u = list(abs(hV-h_u)).index(min(abs(hV-h_u)))   # find a close value in hV
p_u = p_unstableV[idx_u]
plt.text(h_u, p_u-0.01, '$p_u^*$ ', ha='right', va='top', fontsize='x-large')

def make_arrow(x, y, dy):
    plt.annotate('', xy=(x, y), xytext=(x, y+dy), arrowprops=dict(arrowstyle='->'))

if suffix == '_5':
	make_arrow(1-q_0/2,         0.55, -0.1)
	make_arrow((1-q_0+1-q_1)/2,   0.60, -0.1)
	make_arrow((1-q_0+1-q_1)/2,   0.00, +0.07)
	make_arrow((1-q_1+1-q_hat)/2, 0.62, -0.1)
	make_arrow((1-q_1+1-q_hat)/2, 0.10, +0.1)
	make_arrow((1-q_1+1-q_hat)/2, 0.85, +0.1)
	make_arrow((1-q_hat)/2,   0.45, +0.1)

elif suffix == '_7':
	make_arrow(1-q_1/2,         0.55, -0.1)
	make_arrow((1-q_0+1-q_1)/2,   0.51, -0.1)
	make_arrow((1-q_0+1-q_1)/2,   0.90, +0.07)
	make_arrow((1-q_0+1-q_hat)/2, 0.43, -0.1)
	make_arrow((1-q_0+1-q_hat)/2, 0.03, +0.1)
	make_arrow((1-q_0+1-q_hat)/2, 0.75, +0.1)

# decorate the plot

# put critical h points on the x-axis

if h_hat != 0:
    h_locs = [0, h_hat, h_0, h_1, 1]
    h_labs = ['0', r'$\hat{h}$', r'$h_0$', r'$h_1$', '1']
else:
    h_locs = [0, h_0, h_1, 1]
    h_labs = ['0', r'$h_0$', r'$h_1$', '1']
h_zip = sorted(list(zip(h_locs, h_labs)), key = lambda v: v[0])

new_locs, new_labs = zip(*h_zip)
plt.xticks(new_locs, new_labs, fontsize='x-large')

# rest of the decorations

plt.xlabel(r'homophily $h$', fontsize='x-large')
plt.ylabel(r'proportion of Cooperators $p$', fontsize='x-large')
plt.xlim((0, 1))
plt.ylim((plt_ymin, plt_ymax))
plt.legend(loc='upper center', framealpha=1, ncol=2)

plt.tight_layout()
plt.savefig(res_dir + 'isocline_detailed' + suffix + '.pdf')
plt.close()

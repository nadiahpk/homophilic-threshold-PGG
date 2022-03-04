# plot multiple isoclines for different group sizes n

import matplotlib.pyplot as plt
import pandas as pd

# ----------------------

# user parameters
# ---

res_dir = '../../results/leader_driven/'
colours = ['red', 'blue', 'orange', 'black']


# read in the critical points
# ---

# get the critical points
fname = res_dir + 'isoclines_v_n_critpts.csv'
df = pd.read_csv(fname)


# each row of df is one isocline we wish to plot
# ---

height = 3.8 # default 4.8
width = 1.3*height
plt.figure(figsize=(width, height)) # default

plot_lines = list() # for making legends
nV = list()

for row in range(len(df)):

    colour = colours[row]

    # read in the information
    suffix = df.iloc[row]['suffix']
    n = df.iloc[row]['n']
    q_0 = df.iloc[row]['q_0']
    q_1 = df.iloc[row]['q_1']
    q_hat = df.iloc[row]['q_hat']
    p_at_h_hat = df.iloc[row]['p_at_q_hat']

    nV.append(n)

    # read in the corresponding file for the isocline
    fname = res_dir + 'isoclines_v_n' + suffix + '.csv'
    df_iso = pd.read_csv(fname)
    qV = df_iso['q'].values
    p_stableV = df_iso['p_stable'].values
    p_unstableV = df_iso['p_unstable'].values

    # convert to h = 1-q
    hV = 1-qV
    h_0 = 1-q_0
    h_1 = 1-q_1
    h_hat = 1-q_hat

    # plot the isoclines

    l1, = plt.plot(hV, p_stableV,   color=colour, alpha=0.7, lw=3)
    l2, = plt.plot(hV, p_unstableV, color=colour, alpha=0.7, ls='dashed')

    # append for the legend later
    plot_lines.append([l1, l2])

    # plot the p=0 and p=1 continuation of the isoclines

    plt.plot([0, h_0], [0, 0],  color=colour, alpha=0.7, lw=3)
    plt.plot([0, h_1], [1, 1],  color=colour, alpha=0.7, ls='dashed')

    # plot the point where the transition from 2 to 0 interior equilibria occurs (if needed)
    if h_hat != 0:
        plt.scatter([h_hat], [p_at_h_hat], color=colour)




# create legend
# ---

legend1 = plt.legend(plot_lines[3], ['stable', 'unstable'], loc='upper center', framealpha=1, ncol=2)
legend2 = plt.legend([l[0] for l in plot_lines], [r'$n=' + str(n) + r'$' for n in nV], loc='center right')
plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)

# decorate and save the plot
# ---

plt.xlabel(r'homophily $h$', fontsize='x-large')
plt.ylabel(r'proportion of Cooperators $p$', fontsize='x-large')
plt.xlim((0, 1))
plt.ylim((-0.05, 1.15))

plt.tight_layout()
plt.savefig('../../results/leader_driven/isoclines_v_n.pdf')
plt.close()

# make plots showing how the relative abundance in the top x families changes with group size
import pandas as pd
import matplotlib.pyplot as plt


# user parameters
# ---

little_letters = ['(a) leader driven', '(b) members recruit', '(c) members attract']
model_types = [ 'leader_driven', 'members_recruit', 'members_attract' ]
res_dir = '../../results/family_abundance_distn/'

# which ranks to plot
xV = [1, 2, 3]
colours = ['red', 'blue', 'black']
lsV = ['solid', 'dashed', 'dotted']


# make separate plot for each model type
# ---

fig, axs = plt.subplots(1, 3, sharey=True, figsize=(9,4.0))
fig.add_subplot(111, frameon=False)
for panel in range(3):

    model_type = model_types[panel]
    little_letter = little_letters[panel]

    for xi, colour, ls in zip(xV, colours, lsV):

        # read in the data
        fname = res_dir + 'relative_abundances_' + model_type + '_rank_' + str(xi) + '.csv'
        df = pd.read_csv(fname)

        nV = [ int(col_name.split('=')[1]) for col_name in df.columns ]

        # plots
        axs[panel].plot(nV, df.quantile(0.5), color=colour, lw=2, label=str(xi), ls=ls)            # mean
        axs[panel].fill_between(nV, df.quantile(0.05), df.quantile(0.95), color=colour, alpha=0.2) # 90%-ile
        # axs[panel].plot(nV, df.iloc[0], lw=0.5, color=colour) # one example

    if panel == 0:
        axs[panel].legend(loc='best', title='family rank', fontsize='large')
    axs[panel].text(1.0, 1.1, little_letter, fontsize='large')
    axs[panel].set_ylim((0,1))
    axs[panel].set_xlim((2, 600000))
    axs[panel].set_xscale('log')

plt.xticks([0],['hi'], alpha=0) # so there's enough space for the xlabel
plt.yticks([0],['hi i'], alpha=0) # so there's enough space for the ylabel
plt.tight_layout()
plt.xlabel(r'group size $n$', fontsize='x-large')
plt.ylabel(r'proportion of group in family', fontsize='x-large')
plt.savefig(res_dir + 'relative_abundances.pdf')
plt.close()

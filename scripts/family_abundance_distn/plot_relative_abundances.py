# make plots showing how the relative abundance in the top x families changes with group size
import pandas as pd
import matplotlib.pyplot as plt


# user parameters
# ---

model_types = ['members_recruit', 'leader_driven', 'members_attract']
res_dir = '../../results/family_abundance_distn/'

# which ranks to plot
xV = [1, 2, 3]
colours = ['red', 'blue', 'black']
lsV = ['solid', 'dashed', 'dotted']


# make separate plot for each model type
# ---

for model_type in model_types:

    # plot
    height = 3.8 # default 4.8
    width = 1.3*height
    plt.figure(figsize=(width, height)) # default

    for xi, colour, ls in zip(xV, colours, lsV):

        # read in the data
        fname = res_dir + 'relative_abundances_' + model_type + '_rank_' + str(xi) + '.csv'
        df = pd.read_csv(fname)

        nV = [ int(col_name.split('=')[1]) for col_name in df.columns ]

        # plots
        plt.plot(nV, df.quantile(0.5), color=colour, lw=2, label=str(xi), ls=ls)            # mean
        plt.fill_between(nV, df.quantile(0.05), df.quantile(0.95), color=colour, alpha=0.2) # 90%-ile
        # plt.plot(nV, df.iloc[0], lw=0.5, color=colour) # one example


    plt.xscale('log')
    plt.ylim((0,1))
    plt.xlim((2, 600000))
    plt.xlabel(r'group size $n$', fontsize='x-large')
    plt.ylabel(r'proportion of group in family', fontsize='x-large')
    plt.legend(loc='best', title='family rank')
    plt.tight_layout()
    plt.savefig(res_dir + 'relative_abundances_' + model_type + '.pdf')
    plt.close()

# put details in suffix = 5

import pandas as pd
import matplotlib.pyplot as plt

# parameters
# ---

suffix = '_5'
res_dir = '../../results/members_recruit/'
out_dir = '../../results/figs_main/'

# put panel-specific info into dict
little_letters = ['(a) leader driven', '(b) members recruit', '(c) members attract']
res_dirs = [ '../../results/leader_driven/', '../../results/members_recruit/', '../../results/members_attract/' ]

'''
# plot size
height = 3.7 # default 4.8
width = 1.1*height
plt.figure(figsize=(width, height)) # default
'''

fig, axs = plt.subplots(1, 3, sharey=True, figsize=(9,4.0))
fig.add_subplot(111, frameon=False)
for panel in range(3):

    res_dir = res_dirs[panel]
    little_letter = little_letters[panel]

    # read in the data
    # ---

    if panel == 2:

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

    else:

        fname = res_dir + 'critical_v_n' + suffix + '.csv'
        df = pd.read_csv(fname)
        nV = df['nV'].values
        q0V = df['q0V'].values
        q1V = df['q1V'].values
        qhV = df['qhV'].values

        # invert to match the new "homophily" as the variable
        h0V = 1-q0V
        h1V = 1-q1V
        hhV = 1-qhV

    # where cooperators can invade
    axs[panel].plot(nV, h0V, '-o', color='blue')#, label=r'$q_0$')
    axs[panel].fill_between(nV, h0V, 1, alpha=0.25, color='blue', label='cooperators can invade')

    # where defectors can invade
    axs[panel].plot(nV, h1V, '-o', color='red')#, label=r'$q_1$')
    axs[panel].fill_between(nV, h1V, 0, alpha=0.25, color='red', label='defectors can invade')

    # where cooperators cannot persist
    axs[panel].plot(nV, hhV, '-o', color='black')#, label=r'$\hat{q}$')
    hvV = [ min([h0,hh]) for h0, hh in zip(h0V, hhV) ]
    axs[panel].fill_between(nV, hvV, 0, alpha=0.25, color='black', label='cooperators cannot persist')

    # write as text each of the lines' identity

    if h0V[-1] > h1V[-1]:
        h0va = 'bottom'
        h1va = 'top'
    else:
        h0va = 'top'
        h1va = 'bottom'

    axs[panel].annotate(r'  $h_1$ ', (nV[1], h1V[1]), ha='right', va=h1va, fontsize='x-large', color='red')
    axs[panel].annotate(r'  $h_0$', (nV[1], h0V[1]), ha='left', va=h0va, fontsize='x-large', color='blue')
    axs[panel].annotate(r'  $\hat{h}$ ', (nV[2], hhV[2]), ha='right', va='bottom', fontsize='x-large', color='black')

    # decorate with explainers

    # to make the labels to on top of arrows
    bbox = {'fc': 'white', 'edgecolor': 'none'}

    def make_arrow(x, y, dy, color='white'):
        axs[panel].annotate('', xy=(x, y), xytext=(x, y+dy), 
                arrowprops=dict(arrowstyle='|-|', lw=2, shrinkA=0, shrinkB=0, color=color), 
                annotation_clip=False)

    offset = 0
    posn = 5
    if panel == 0:

        make_arrow(nV[posn]-offset, h1V[posn], -h1V[posn]), 
        axs[panel].text(nV[posn]-offset, h1V[posn]/2, ' D can invade \n & persist ', ha='center', va='center', bbox=bbox)
        make_arrow(nV[posn]-offset, h1V[posn], 1-h1V[posn]), 
        axs[panel].text(nV[posn]-offset, h1V[posn] + (1-h1V[posn])/2, ' D cannot \n invade nor \n persist ', ha='center', va='center', bbox=bbox)

    elif panel == 2:

        #posn = 5
        make_arrow(nV[posn]+offset, h0V[posn], 1-h0V[posn])
        axs[panel].text(nV[posn]+offset, h0V[posn] + (1-h0V[posn])/2, ' C can \n invade ', ha='center', va='center', bbox=bbox)
        make_arrow(nV[posn]+offset, h0V[posn], -h0V[posn])
        axs[panel].text(nV[posn]+offset, h0V[posn]/2, ' C cannot \n invade ', ha='center', va='center', bbox=bbox)

    else: # panel == 1:

        #posn = 6
        make_arrow(nV[posn], hhV[posn], 1-hhV[posn])
        axs[panel].text(nV[posn], hhV[posn] + (1-hhV[posn])/2, ' C can \n persist ', ha='center', va='center', bbox=bbox)

        make_arrow(nV[posn], 0, hhV[posn])
        axs[panel].text(nV[posn], hhV[posn]/2, ' C cannot \n persist ', ha='center', va='center', bbox=bbox)

    axs[panel].text(3.0, 1.1, little_letter, fontsize='large')
    axs[panel].set_xticks([4, 18])
    axs[panel].set_xlim((nV[0]-0.5, nV[-1]+0.5))
    axs[panel].set_ylim((-0.02, 1.02))

# finish off

plt.xticks([0],['hi'], alpha=0) # so there's enough space for the xlabel
plt.yticks([0],['hi i'], alpha=0) # so there's enough space for the ylabel
plt.ylabel(r'homophily $h$', fontsize='x-large')
plt.xlabel(r'group size $n$', fontsize='x-large')
plt.tight_layout()
plt.savefig(out_dir + 'critical_v_n' + suffix + '.pdf')
plt.close()

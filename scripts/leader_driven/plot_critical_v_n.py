# put details in suffix = 5

import pandas as pd
import matplotlib.pyplot as plt

# parameters
# ---

#suffix = '_6'
suffix = '_7'
suffix = '_5'
#suffix = '_4242'
res_dir = '../../results/leader_driven/'


# read in the data
# ---

fname = res_dir + 'critical_v_n' + suffix + '.csv'
df = pd.read_csv(fname)
nV = df['nV'].values
q0V = df['q0V'].values
q1V = df['q1V'].values
qhV = df['qhV'].values

# invert to match the new "homophily" as the variable
q0V = 1-q0V
q1V = 1-q1V
qhV = 1-qhV

# plot size
height = 3.7 # default 4.8
width = 1.1*height
plt.figure(figsize=(width, height)) # default

# where cooperators can invade
plt.plot(nV, q0V, '-o', color='blue')#, label=r'$q_0$')
plt.fill_between(nV, q0V, 1, alpha=0.25, color='blue', label='cooperators can invade')

# where defectors can invade
plt.plot(nV, q1V, '-o', color='red')#, label=r'$q_1$')
plt.fill_between(nV, q1V, 0, alpha=0.25, color='red', label='defectors can invade')

# where cooperators cannot persist
plt.plot(nV, qhV, '-o', color='black')#, label=r'$\hat{q}$')
qvV = [ min([q0,qh]) for q0, qh in zip(q0V, qhV) ]
plt.fill_between(nV, qvV, 0, alpha=0.25, color='black', label='cooperators cannot persist')

# write as text each of the lines' identity

if q0V[-1] > q1V[-1]:
    q0va = 'bottom'
    q1va = 'top'
else:
    q0va = 'top'
    q1va = 'bottom'

plt.annotate(r'  $h_1$ ', (nV[1], q1V[1]), ha='right', va=q1va, fontsize='x-large', color='red')
plt.annotate(r'  $h_0$', (nV[1], q0V[1]), ha='left', va=q0va, fontsize='x-large', color='blue')
plt.annotate(r'  $\hat{h}$ ', (nV[2], qhV[2]), ha='right', va='bottom', fontsize='x-large', color='black')

if False:

    # decorate with explainers

    # to make the labels to on top of arrows
    bbox = {'fc': 'white', 'edgecolor': 'none', 'pad': 2}

    def make_arrow(x, y, dy, color='white'):
        plt.annotate('', xy=(x, y), xytext=(x, y+dy), 
                arrowprops=dict(arrowstyle='|-|', lw=2, shrinkA=0, shrinkB=0, color=color), 
                annotation_clip=False)

    offset = 0
    posn = 4
    make_arrow(nV[posn]-offset, q1V[posn], -q1V[posn]), 
    plt.text(nV[posn]-offset, q1V[posn]/2, ' D can invade & persist ', ha='center', va='center', bbox=bbox, rotation=90)
    make_arrow(nV[posn]-offset, q1V[posn], 1-q1V[posn]), 
    plt.text(nV[posn]-offset, q1V[posn] + (1-q1V[posn])/2, ' D cannot \n invade nor \n persist ', ha='center', va='center', bbox=bbox, rotation=90)

    posn = 5
    make_arrow(nV[posn]+offset, q0V[posn], 1-q0V[posn])
    plt.text(nV[posn]+offset, q0V[posn] + (1-q0V[posn])/2, ' C can invade ', ha='center', va='center', bbox=bbox, rotation=90)
    make_arrow(nV[posn]+offset, q0V[posn], -q0V[posn])
    plt.text(nV[posn]+offset, q0V[posn]/2, ' C cannot invade ', ha='center', va='center', bbox=bbox, rotation=90)

    posn = 6
    make_arrow(nV[posn], qhV[posn], 1-qhV[posn])
    plt.text(nV[posn], qhV[posn] + (1-qhV[posn])/2, ' C can persist ', ha='center', va='center', bbox=bbox, rotation=90)

    make_arrow(nV[posn], 0, qhV[posn])
    plt.text(nV[posn], qhV[posn]/2, ' C cannot \n persist ', ha='center', va='center', bbox=bbox, rotation=90)


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

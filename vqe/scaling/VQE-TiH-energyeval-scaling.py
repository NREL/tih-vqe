#!/usr/bin/env python


import os,sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

dpi = 300

#initialize plot
matplotlib.rc('font',size=14)
matplotlib.rc('axes',titlesize = 14)
matplotlib.rc('figure',titlesize = 14)
matplotlib.rc('xtick',labelsize=14)
matplotlib.rc('ytick',labelsize=14)
matplotlib.rc('legend',fontsize = 14)
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)

#initialize plot, set number of subplots
fig = plt.figure(figsize=(6,5),dpi=dpi)
gs1 = gridspec.GridSpec(1,1)
axmain = fig.add_subplot(gs1[0],frameon=False)
gs2 = gridspec.GridSpec(1,1,height_ratios=[1],width_ratios=[1])

ax = []  #list of ax
for i in range(0,1):
    ax.append(fig.add_subplot(gs2[i]))
gs2.tight_layout(fig,rect=[0.06,0.02,1,1])
gs2.update(hspace=0.1)
gs2.update(wspace=0.1)
gs1.tight_layout(fig,rect=[0,0,1,0.96])

axmain.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
axmain.set_xticklabels([])
axmain.set_yticklabels([])
axmain.grid(False)
axmain.set_zorder(1)
axmain.set_visible(False)
#axmain.yaxis.labelpad = 5

#data
#processors
p = np.array([2,4,8,16,36,72,144,216])
#time per step (seconds)
t = np.array([2134.90,1108.57,564.39,287.69,126.17,65.79,34.68,26.58])
Nfirsttoplot = 5   #only plot first many of these values in second plot

timepereval = 1/t*3600.0   #convert to evaluations per hour

#actually plot data
xmin = 0
xmax = 220
xs = 50
ymin = 0
ymax = 150
ys = 30

i = 0  #first subplot
ax[i].set_zorder(2)
line1, = ax[i].plot(p,timepereval, color='dimgray',markeredgecolor='k', zorder = 6,
                    linewidth=2,linestyle='-',marker='o',markersize=10,clip_on=False)
ax[i].axhline(y=0, color='k', linestyle='-',zorder=5)
ax[i].axvline(x=0, color='k', linestyle='-',zorder=5)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_xticklabels(range(xmin,xmax+1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+ys/2,ys))
ax[i].set_yticklabels(range(ymin,ymax+ys,ys))
ax[i].set_ylabel('TiH energy evaluations per hour')
ax[i].set_xlabel('Number of parallel processors')
#ax[i].legend((line1),('Parity mapping'),edgecolor='gainsboro',loc='upper left',fontsize=12) #make legend, need at least 2 lines
#can use custom lines if need to
#custom_lines = [Line2D([0], [0], color='gray', lw=3),
#                Line2D([0], [0], color='gray', lw=3, ls='--')]
#leg = plt.figlegend(custom_lines,['DFT','RPA@DFT'],loc='upper center',ncol=2,frameon=False,bbox_to_anchor=[0.5,0.98])

for i in range(0,1):
    [ax[i].spines[x].set_linewidth(3) for x in ax[i].spines]
#    [x.set_zorder(3) for x in ax[i].spines.itervalues()]
    ax[i].tick_params(direction = 'out',width = 3,length=5)
    ax[i].tick_params(top='off',right='off')
    ax[i].grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax[i].set_axisbelow(True)

plt.tight_layout()
plt.savefig('VQE-TiH-energyeval-scaling.png')
plt.show()



p = p[:Nfirsttoplot]
timepereval = timepereval[:Nfirsttoplot]

#initialize plot, set number of subplots
fig = plt.figure(figsize=(6,5),dpi=dpi)
gs1 = gridspec.GridSpec(1,1)
axmain = fig.add_subplot(gs1[0],frameon=False)
gs2 = gridspec.GridSpec(1,1,height_ratios=[1],width_ratios=[1])

ax = []  #list of ax
for i in range(0,1):
    ax.append(fig.add_subplot(gs2[i]))
gs2.tight_layout(fig,rect=[0.06,0.02,1,1])
gs2.update(hspace=0.1)
gs2.update(wspace=0.1)
gs1.tight_layout(fig,rect=[0,0,1,0.96])

axmain.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
axmain.set_xticklabels([])
axmain.set_yticklabels([])
axmain.grid(False)
axmain.set_zorder(1)
axmain.set_visible(False)
#axmain.yaxis.labelpad = 5

#actually plot data
xmin = 0
xmax = 50
xs = 10
ymin = 0
ymax = 30
ys = 5

i = 0  #first subplot
ax[i].set_zorder(2)
line1, = ax[i].plot(p,timepereval, color='k',markeredgecolor='k', zorder = 6,
                    linewidth=2,linestyle='-',marker='o',markersize=10,clip_on=False)
ax[i].axhline(y=0, color='k', linestyle='-',zorder=5)
ax[i].axvline(x=0, color='k', linestyle='-',zorder=5)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_xticklabels(range(xmin,xmax+1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+ys/2,ys))
ax[i].set_yticklabels(range(ymin,ymax+ys,ys))
ax[i].set_ylabel('TiH energy evaluations per hour')
ax[i].set_xlabel('Number of parallel processors')
#ax[i].legend((line1),('Parity mapping'),edgecolor='gainsboro',loc='upper left',fontsize=12) #make legend, need at least 2 lines
#can use custom lines if need to
#custom_lines = [Line2D([0], [0], color='gray', lw=3),
#                Line2D([0], [0], color='gray', lw=3, ls='--')]
#leg = plt.figlegend(custom_lines,['DFT','RPA@DFT'],loc='upper center',ncol=2,frameon=False,bbox_to_anchor=[0.5,0.98])

for i in range(0,1):
    [ax[i].spines[x].set_linewidth(3) for x in ax[i].spines]
#    [x.set_zorder(3) for x in ax[i].spines.itervalues()]
    ax[i].tick_params(direction = 'out',width = 3,length=5)
    ax[i].tick_params(top='off',right='off')
    ax[i].grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax[i].set_axisbelow(True)

plt.tight_layout()
plt.savefig('VQE-TiH-energyeval-scaling-small.png')
plt.show()


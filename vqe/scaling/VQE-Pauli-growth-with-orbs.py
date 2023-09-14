#!/usr/bin/env python


import os,sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

#initialize plot
matplotlib.rc('font',size=16)
matplotlib.rc('axes',titlesize = 16)
matplotlib.rc('figure',titlesize = 16)
matplotlib.rc('xtick',labelsize=16)
matplotlib.rc('ytick',labelsize=16)
matplotlib.rc('legend',fontsize = 16)
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)

dpi = 300

#initialize plot, set number of subplots
fig = plt.figure(figsize=(6,5),dpi=dpi)
gs1 = gridspec.GridSpec(1,1)
axmain = fig.add_subplot(gs1[0],frameon=False)
gs2 = gridspec.GridSpec(1,1,height_ratios=[1],width_ratios=[1])

ax = []  #list of ax
for i in range(0,1):
    ax.append(fig.add_subplot(gs2[i]))
gs2.tight_layout(fig,rect=[0.06,0.05,1,1])
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
x = [4,6,12,16]  #spin orbitals
orig = [5,100,631,2913]
greedy = [2,25,175,611]
x = [4,6,10,12,16]
orig = [5,100,276,631,2329]
greedy = [2,25,76,175,566]

#actually plot data
xmin = 0
xmax = 20
xs = 5
ymin = 0
ymax = 3000
ys = 500

i = 0  #first subplot
ax[i].set_zorder(2)
line1, = ax[i].plot(x,orig, color='k', zorder = 6,linewidth=2,linestyle='--',marker='o',markersize=10,clip_on=False)
line2, = ax[i].plot(x,greedy, color='xkcd:emerald', marker='o',linestyle='--',zorder = 6,linewidth=2,markersize=10,clip_on=False)
ax[i].axhline(y=0, color='k', linestyle='-',zorder=5)
ax[i].axvline(x=0, color='k', linestyle='-',zorder=5)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_xticklabels(range(xmin,xmax+1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+0.1,ys))
ax[i].set_yticklabels(range(ymin,ymax+1,ys))
ax[i].set_ylabel('Number of Pauli strings to measure')
ax[i].set_xlabel('Number of spin orbitals')
ax[i].legend((line1,line2),('Parity mapping','Parity mapping + greedy'),edgecolor='gainsboro',loc='upper left',fontsize=14) #make legend, need at least 2 lines
#can use custom lines if need to
#custom_lines = [Line2D([0], [0], color='gray', lw=3),
#                Line2D([0], [0], color='gray', lw=3, ls='--')]
#leg = plt.figlegend(custom_lines,['DFT','RPA@DFT'],loc='upper center',ncol=2,frameon=False,bbox_to_anchor=[0.5,0.98])
ax[i].annotate('H$_2$ 2 qubits',xytext=(x[0],550),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('(Li,Na,K)H 4 qubits',xytext=(x[1],1000),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('(Li,Na,K)H 8 qubits',xytext=(x[2],1200),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('TiH 10 qubits',xytext=(x[3],1350),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('TiH 14 qubits',xytext=(x[4],1600),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')

for i in range(0,1):
    [ax[i].spines[x].set_linewidth(3) for x in ax[i].spines]
#    [x.set_zorder(3) for x in ax[i].spines.itervalues()]
    ax[i].tick_params(direction = 'out',width = 3,length=5)
    ax[i].tick_params(top='off',right='off')
    ax[i].grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax[i].set_axisbelow(True)

plt.tight_layout()
plt.savefig('VQE-Pauli-growth-with-orbs.png')
plt.show()


#initialize plot, set number of subplots
fig = plt.figure(figsize=(6,5),dpi=dpi)
gs1 = gridspec.GridSpec(1,1)
axmain = fig.add_subplot(gs1[0],frameon=False)
gs2 = gridspec.GridSpec(1,1,height_ratios=[1],width_ratios=[1])

ax = []  #list of ax
for i in range(0,1):
    ax.append(fig.add_subplot(gs2[i]))
gs2.tight_layout(fig,rect=[0.06,0.05,1,1])
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
xmax = 20
xs = 5
ymin = 0
ymax = 3000
ys = 500

i = 0  #first subplot
ax[i].set_zorder(2)
line1, = ax[i].plot(x,orig, color='k', zorder = 6,linewidth=2,linestyle='--',marker='o',markersize=10,clip_on=False)
#line2, = ax[i].plot(x,greedy, color='xkcd:emerald', marker='o',linestyle='--',zorder = 6,linewidth=2,markersize=10,clip_on=False)
ax[i].axhline(y=0, color='k', linestyle='-',zorder=5)
ax[i].axvline(x=0, color='k', linestyle='-',zorder=5)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_xticklabels(range(xmin,xmax+1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+0.1,ys))
ax[i].set_yticklabels(range(ymin,ymax+1,ys))
ax[i].set_ylabel('Number of Pauli strings to measure')
ax[i].set_xlabel('Number of spin orbitals')
#ax[i].legend((line1,line2),('Parity mapping','Parity mapping + greedy'),edgecolor='gainsboro',loc='upper left',fontsize=14) #make legend, need at least 2 lines
ax[i].legend(tuple([line1]),tuple(['Parity mapping']),edgecolor='gainsboro',loc='upper left',fontsize=14) #make legend, need at least 2 lines

#can use custom lines if need to
#custom_lines = [Line2D([0], [0], color='gray', lw=3),
#                Line2D([0], [0], color='gray', lw=3, ls='--')]
#leg = plt.figlegend(custom_lines,['DFT','RPA@DFT'],loc='upper center',ncol=2,frameon=False,bbox_to_anchor=[0.5,0.98])
ax[i].annotate('H$_2$ 2 qubits',xytext=(x[0],550),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('(Li,Na,K)H 4 qubits',xytext=(x[1],1000),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('(Li,Na,K)H 8 qubits',xytext=(x[2],1200),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('TiH 10 qubits',xytext=(x[3],1350),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')
ax[i].annotate('TiH 14 qubits',xytext=(x[4],1600),xy=(0,0),annotation_clip=False,fontsize=14,rotation=90,ha='center',va='center')

for i in range(0,1):
    [ax[i].spines[x].set_linewidth(3) for x in ax[i].spines]
#    [x.set_zorder(3) for x in ax[i].spines.itervalues()]
    ax[i].tick_params(direction = 'out',width = 3,length=5)
    ax[i].tick_params(top='off',right='off')
    ax[i].grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax[i].set_axisbelow(True)

plt.tight_layout()
plt.savefig('VQE-Pauli-growth-with-orbs-2.png')
plt.show()


#!/usr/bin/env python

import os,sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

#initialize plot
matplotlib.rc('font',size=14)
matplotlib.rc('axes',titlesize = 14)
matplotlib.rc('figure',titlesize = 14)
matplotlib.rc('xtick',labelsize=12)
matplotlib.rc('ytick',labelsize=12)
matplotlib.rc('legend',fontsize = 12)
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)

#initialize plot, set number of subplots
fig = plt.figure(figsize=(6,5),dpi=300)
gs1 = gridspec.GridSpec(1,1)
axmain = fig.add_subplot(gs1[0],frameon=False)
gs2 = gridspec.GridSpec(2,2,height_ratios=[1,1],width_ratios=[1,1])

ax = []  #list of ax
for i in range(0,4):
    ax.append(fig.add_subplot(gs2[i]))
gs2.tight_layout(fig,rect=[0.04,0.1,1,0.94])
gs2.update(hspace=0.3)
gs2.update(wspace=0.45)
gs1.tight_layout(fig,rect=[0,0,1,0.96])

axmain.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
axmain.set_xticklabels([])
axmain.set_yticklabels([])
axmain.grid(False)
axmain.set_zorder(1)
#axmain.yaxis.labelpad = 5
axmain.set_visible(False)

xmin = 0
xmax = 4
xs = 1

exactc = 'k'
uccsc = 'dodgerblue'
uccsdc = 'orange'

#data
#xlsx not supported for some reason??? old version?
filename = 'dissociation-curves.xls'
x = np.arange(0.6,3.05,0.1)

mysheet = 'LiH 3'
i = 0  #first subplot
data = pd.read_excel(filename,sheet_name=mysheet,skiprows=1)
data = data.iloc[:25,:]
exact = data['Exact']
uccs = data['UCCS']
uccsd = data['UCCSD']
ymin = -8.0
ymax = -7.2
ymin = -0.2
ymax = 0.7
ys = 0.2
ax[i].set_zorder(2)
line1, = ax[i].plot(x,exact-exact.iloc[-1], color=exactc, zorder = 8,linewidth=1)
line2, = ax[i].plot(x,uccs-exact.iloc[-1], color=uccsc, marker='o',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
line3, = ax[i].plot(x,uccsd-exact.iloc[-1], color=uccsdc, marker='s',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+0.1,ys))
ax[i].set_yticklabels([round(x,2) for x in np.arange(ymin,ymax+1,ys)])
ax[i].set_title('LiH')

mysheet = 'NaH 4'
i = 1  #second subplot
data = pd.read_excel(filename,sheet_name=mysheet,skiprows=1)
data = data.iloc[:25,:]  #don't plot above 3 Ang
exact = data['Exact']
uccs = data['UCCS']
uccsd = data['UCCSD']
ymin = -160.4
ymax = -159.6
ymin = -0.2
ymax = 0.7
ys = 0.2
ax[i].set_zorder(2)
line1, = ax[i].plot(x,exact-exact.iloc[-1], color=exactc, zorder = 8,linewidth=1)
line2, = ax[i].plot(x,uccs-exact.iloc[-1], color=uccsc, marker='o',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
line3, = ax[i].plot(x,uccsd-exact.iloc[-1], color=uccsdc, marker='s',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+0.1,ys))
ax[i].set_yticklabels([round(x,2) for x in np.arange(ymin,ymax+1,ys)])
ax[i].set_title('NaH')

mysheet = 'KH 5'
i = 2  #third subplot
data = pd.read_excel(filename,sheet_name=mysheet,skiprows=1)
data = data.iloc[:25,:]
exact = data['Exact']
uccs = data['UCCS']
uccsd = data['UCCSD']
ymin = -593.7
ymax = -592.9
ymin = -0.2
ymax = 0.7
ys = 0.2
ax[i].set_zorder(2)
line1, = ax[i].plot(x,exact-exact.iloc[-1], color=exactc, zorder = 8,linewidth=1)
line2, = ax[i].plot(x,uccs-exact.iloc[-1], color=uccsc, marker='o',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
line3, = ax[i].plot(x,uccsd-exact.iloc[-1], color=uccsdc, marker='s',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+0.1,ys))
ax[i].set_yticklabels([round(x,2) for x in np.arange(ymin,ymax+1,ys)])
ax[i].set_title('KH')

mysheet = 'TiH-combinedspins-lg'
i = 3  #fourth subplot
data = pd.read_excel(filename,sheet_name=mysheet,skiprows=1)
data = data.iloc[:26,:]
exact = data['Exact']
uccs = data['UCCS']
uccsd = data['UCCSD']
ymin = -840.1
ymax = -839.3
ymin = -0.2
ymax = 0.7
ys = 0.2
ax[i].set_zorder(2)
line1, = ax[i].plot(x,exact-exact.iloc[-1], color=exactc, zorder = 8,linewidth=1)
line2, = ax[i].plot(x,uccs-exact.iloc[-1], color=uccsc, marker='o',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
line3, = ax[i].plot(x,uccsd-exact.iloc[-1], color=uccsdc, marker='s',linestyle='-',zorder = 6,linewidth=1,markersize=3,alpha=1)
ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[i].set_xticks(np.arange(xmin,xmax+0.1,xs))
ax[i].set_yticks(np.arange(ymin,ymax+0.1,ys))
ax[i].set_yticklabels([round(x,2) for x in np.arange(ymin,ymax+1,ys)])
ax[i].set_title('TiH')

for i in range(0,4):
    #ax[i].axhline(y=0, color='k', linestyle='-',zorder=5)
    #ax[i].axvline(x=0, color='k', linestyle='-',zorder=5)
    if i == 2 or i == 3:
        ax[i].set_xticklabels(np.arange(xmin,xmax+0.1,xs))
        ax[i].set_xlabel('Bond length ($\AA$)')
    else:
        ax[i].set_xticklabels([])
    if i == 0 or i == 2:
        ax[i].set_ylabel('E - E$_\mathrm{dissc}$ (Ha)')
    ax[i].legend((line1,line2,line3),('Exact','UCCS','UCCSD'),edgecolor='gainsboro',loc='upper right',fontsize=10) #make legend, need at least 2 lines

for i in range(0,4):
    [ax[i].spines[x].set_linewidth(2) for x in ax[i].spines]
    [ax[i].spines[x].set_zorder(30) for x in ax[i].spines]
#    [ax[i].spines[x].set_zorder(3) for x in ax[i].spines]
    ax[i].tick_params(direction = 'out',width = 2,length=4)
    #ax[i].tick_params(top='off',right='off')
    ax[i].grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax[i].set_axisbelow(True)

plt.tight_layout()
plt.savefig('bond-dissociation-hydrides.png')
plt.show()


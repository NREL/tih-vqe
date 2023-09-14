#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os,sys
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import json
import math
import re
import matplotlib.colors as myc


matplotlib.rc('font',size=16)
matplotlib.rc('axes',titlesize = 16)
matplotlib.rc('figure',titlesize = 16)
matplotlib.rc('xtick',labelsize=14)
matplotlib.rc('ytick',labelsize=14)
matplotlib.rc('legend',fontsize = 14)
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)

def plotdata(ax,xdata,ydata,title,i,yrng):
    mincol = np.argmin(np.min(ydata,axis=0))
    if 'LiH' in title:
        lines = []
        for c in range(2):
            mfc = myc.to_rgba(colors[c], 0)
            if c == mincol:
                mfc = colors[c]
            line1, = ax.plot(xdata,ydata[:,c]-ydata[-1,1],markeredgecolor=colors[c],markerfacecolor=mfc,color=colors[c],marker='o',markersize=5,zorder=5,clip_on=False)
            lines.append(line1)
        lines = tuple(lines)
        labels = tuple(['S = 1','S = 3'])
    elif 'TiH' in title:
        lines = []
        for c in range(3):
            mfc = myc.to_rgba(colors[c], 0)
            if c == mincol:
                mfc = colors[c]
            line1, = ax.plot(xdata,ydata[:,c]-ydata[-1,1],markeredgecolor=colors[c],markerfacecolor=mfc,color=colors[c],marker='o',markersize=5,zorder=5,clip_on=False)
            lines.append(line1)
        lines = tuple(lines)
        labels = tuple(['S = 2','S = 4','S = 6'])

    ax.legend(lines,labels)
    ax.axhline(y=0,color='k',linewidth=0.5,zorder=1)

    ax.set_xticks(np.arange(1,2.6,0.5))
    ax.set_xlim([1,2.5])
    if i in [2,3]:
        ax.set_xlabel('Bond length ($\AA$)')
    if i in [0,2]:
        ax.set_ylabel('E - E$_\mathrm{dissc}$ (Ha)')
    ax.set_title(title,y=1.02)

    ax.set_ylim(yrng)

    [ax.spines[x].set_linewidth(3) for x in ax.spines]
    [ax.spines[x].set_zorder(3) for x in ax.spines]
    ax.tick_params(direction = 'out',width = 3,length=5)
    ax.grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax.set_axisbelow(True)

    return ax

filename = 'dissociation-curves-gaussian.xls'

data = pd.read_excel(filename,sheet_name='lih',header=None)

bondlengths = np.arange(1,2.51,0.1)

lihsto3g = data.iloc[7:23,1:3].to_numpy()
lihaugccq = data.iloc[7:23,17:19].to_numpy()

colors = ['m','dodgerblue','orange']

fig = plt.figure(figsize=(7,7),dpi=300)
i = 0
ax = fig.add_subplot(2,2,i+1)
ax = plotdata(ax,bondlengths,lihsto3g,'LiH\nLi basis: STO-3G',i,[-0.2,0.2])

i = 1
ax = fig.add_subplot(2,2,i+1)
ax = plotdata(ax,bondlengths,lihaugccq,'LiH\nLi basis: aug-cc-pVQZ',i,[-0.2,0.2])


data = pd.read_excel(filename,sheet_name='tih-ti_only_basis',header=None)
tihsto3g = data.iloc[7:23,1:4].to_numpy()
tihaugccq = data.iloc[7:23,17:20].to_numpy()
i = 2
ax = fig.add_subplot(2,2,i+1)
ax = plotdata(ax,bondlengths,tihsto3g,'TiH\nTi basis: STO-3G',i,[-0.2,0.55])

i = 3
ax = fig.add_subplot(2,2,i+1)
ax = plotdata(ax,bondlengths,tihaugccq,'TiH\nTi basis: aug-cc-pVQZ',i,[-0.2,0.55])



plt.tight_layout()
plt.savefig('gaussian-curves.png')
plt.show()



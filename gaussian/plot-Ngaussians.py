#!/usr/bin/env python3
# -*- coding: utf-8 -*-



import os,sys
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
import pandas
import json
import math
import re


matplotlib.rc('font',size=16)
matplotlib.rc('axes',titlesize = 16)
matplotlib.rc('figure',titlesize = 16)
matplotlib.rc('xtick',labelsize=14)
matplotlib.rc('ytick',labelsize=14)
matplotlib.rc('legend',fontsize = 14)
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)


def readfile(filename):
    #read a file into a list of strings
    f = open(filename,'r')
    tempfile = f.readlines()
    f.close()
    return tempfile

def writefile(filename,tempfile):
    #write tempfile (list of strings) to filename
    f = open(filename,'w')
    f.writelines(tempfile)
    f.close()

def readjson(filename):
    #read json in dict
    with open(filename) as json_file:
        data = json.load(json_file)
    return data

def writejson(filename,dictdata):
    #write dict data to json
    with open(filename,'w') as outfile:
        json.dump(dictdata, outfile)

def find_key(key_input, tempfile):
    #finds line where key occurs in stored input, last instance
    key_input = str(key_input)
    line = len(tempfile)                  #default to end
    for i in range(0,len(tempfile)):
        if key_input in tempfile[i]:
            line = i
    return line

def delta_func(E, a, mu = 0):
    delta = [(1/abs(a)/(2*3.14159)**0.5)*math.exp(-0.5*((x-mu)/a)**2) for x in E]
    return np.array(delta)


mats = ['lih','nah','kh','tih']
files = ['ccsd-sto3g','ccsd-321g','ccsd-631g','ccsd-ccd','ccsd-augccd','ccsd-ccq','ccsd-augccq']
#files = ['ccsd-sto3g','ccsd-631g','ccsd-aug']
#files = ['ccsd-sto3g','ccsd-631g','ccsd-aug','b3lyp-sto3g','b3lyp-631g','b3lyp-aug']
mults = {'tih':4,'nah':1,'kh':1,'lih':1}

Ngauss = {m:{b:np.nan for b in files} for m in mats}
for m in mats:
    for f in files:
        outfile = readfile('experimental-length/'+m+'-'+str(mults[m])+'-'+f)
        try:
            line = find_key('basis functions,',outfile)
            Ng = int(outfile[line].split()[3])
            Ngauss[m][f] = Ng
        except:
            print(m,f)


mattitles = {'tih':'TiH','kh':'KH','nah':'NaH','lih':'LiH'}
basistitles = {'ccsd-sto3g':['STO-3G'],'ccsd-631g':['6-31G'],'ccsd-321g':['3-21G'],
             'ccsd-6311g':['6-311G'],'ccsd-631gdp':['6-31(d\',p\')'],
             'ccsd-ccd':['cc-pVDZ'],'ccsd-augccd':['aug-cc-pVDZ'],
             'ccsd-ccq':['cc-pVQZ'],'ccsd-augccq':['aug-cc-pVQZ'],
          'b3lyp-631g':['6-31G'],'b3lyp-sto3g':['STO-3G'],'b3lyp-aug':['aug-cc-pVQZ']}
xlabels = [basistitles[b][0] for b in files]

xmin = 0
xmax = 8
xs = 1
xrng = np.arange(xmin,xmax+xs/2,xs)
ymin = 0
ymax = 600

fig = plt.figure(figsize=(4,4),dpi=200)
ax = fig.add_subplot(111)

colors = {'lih':'m','nah':'dodgerblue','kh':'xkcd:emerald','tih':'orange'}
lines = []
labels = [mattitles[m] for m in mats]
for m in mats:
    xvals = range(1,len(files)+1)
    yvals = [Ngauss[m][b] for b in files]

    line1, = ax.plot(xvals,yvals,marker='o',linewidth=2,zorder=5,alpha=1,color=colors[m])
    lines.append(line1)

ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_xticks(range(1,xmax))
ax.set_xticklabels(xlabels,rotation=90)
ax.set_xlabel('Cation basis set')
ax.set_ylabel('Number primitive\nGaussian functions')
ax.legend(tuple(lines),tuple(labels),loc='upper left',ncol=2,fontsize=10,framealpha=0.95 )

[ax.spines[x].set_linewidth(3) for x in ax.spines]
[ax.spines[x].set_zorder(30) for x in ax.spines]
ax.tick_params(direction = 'out',width = 3,length=5)
ax.grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig('Ngaussians.png')


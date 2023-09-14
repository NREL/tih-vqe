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
import operator


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
#mats = ['kh']
matlabels = ['LiH','NaH','KH','TiH']
mults = {'lih':[1],'nah':[1],'kh':[1],'tih':[2,4]}
mults = {'lih':[1],'nah':[1],'kh':[1],'tih':[4]}
basissets = ['sto3g','321g','631g','ccd','augccd','ccq','augccq']
basissetlabels = ['STO-3G','3-21G','6-31G','cc-pVDZ','aug-cc-pVDZ','cc-pVQZ','aug-cc-pVQZ']

ev2ryd = 13.605698066
ryd2ha = 2.0

fig = plt.figure(figsize=(4,5),dpi=200)
ax = fig.add_subplot(111)
allbasis = {x:[] for x in mats}
allgaps = {x:[] for x in mats}
alldata = {x:{'gap':[],'HOMO':[],'LUMO':[],'iHOMO':[],'iLUMO':[],
              'eigs':{'up':[],'dn':[]},'Nspin':1,'fermi':[],
              'error':[],'basis':[],'Nbasis':[]} for x in mats}
colors = ['m','dodgerblue','xkcd:emerald','orange']
lines = []
for imat,mat in enumerate(mats):
    for m in mults[mat]:
        for ib,b in enumerate(basissets):
            data = readfile('experimental-length/'+mat+'-'+str(m)+'-ccsd-'+b)

            alldata[mat]['basis'].append(b)
            line = find_key('Error termination',data)
            if line != len(data):
                eg = np.nan
                allbasis[mat].append(np.nan)
                allgaps[mat].append(np.nan)
                alldata[mat]['error'].append(True)

                alldata[mat]['gap'].append(np.nan)
                alldata[mat]['Nbasis'].append(np.nan)
                alldata[mat]['HOMO'].append(np.nan)
                alldata[mat]['LUMO'].append(np.nan)
                alldata[mat]['iHOMO'].append(np.nan)
                alldata[mat]['iLUMO'].append(np.nan)
                alldata[mat]['fermi'].append(np.nan)
                alldata[mat]['eigs']['up'].append([np.nan])
                if alldata[mat]['Nspin'] == 2:
                    alldata[mat]['eigs']['dn'].append([np.nan])
                print('error in out file',mat,m,b)
                continue
            alldata[mat]['error'].append(False)

            line = find_key('NBasis',data)
            Nbasis = data[line].split()[1]
            Nbasis = int(Nbasis)
            alldata[mat]['Nbasis'].append(Nbasis)

            HOMO = -100000
            LUMO = 100000
            iHOMO = [None,None]
            iLUMO = [None,None]
            startline = find_key('Orbital symmetries',data)
            while 'Alpha  occ.' not in data[startline]:
                startline += 1  #go down to eig info
            endline = find_key('Condensed to atoms',data)
            firstalpha = True
            firstbeta = True
            tempeigsup = []
            tempeigsdn = []
            for line in range(startline,endline):
                temp = data[line].split()
                if temp[0] == 'Alpha' and firstalpha:
                    ieig = 0
                    isp = 0
                    firstalpha = False
                elif temp[0] == 'Beta' and firstbeta:
                    ieig = 0
                    isp = 1
                    alldata[mat]['Nspin'] = 2
                    firstbeta = False
                if temp[1] == 'occ.':
                    temp2 = data[line].split('--')[1]
                    eigs = [float(y)*ev2ryd*ryd2ha for y in temp2.split()]
                    if isp == 0:
                        tempeigsup.extend(eigs)
                    else:
                        tempeigsdn.extend(eigs)
                    imaxe, maxe = max(enumerate(eigs), key=operator.itemgetter(1))
                    if maxe > HOMO:
                        HOMO = maxe
                        iHOMO = [isp,ieig+imaxe]
                elif temp[1] == 'virt.':
                    temp2 = data[line].split('--')[1]
                    eigs = [float(y)*ev2ryd*ryd2ha for y in temp2.split()]
                    if isp == 0:
                        tempeigsup.extend(eigs)
                    else:
                        tempeigsdn.extend(eigs)
                    imine, mine = min(enumerate(eigs), key=operator.itemgetter(1))
                    if mine < LUMO:
                        LUMO = mine
                        iLUMO = [isp,ieig+imine]
                ieig += len(eigs)
                #print(ieig)
            eg = LUMO - HOMO
            allbasis[mat].append(Nbasis)
            allgaps[mat].append(eg)
            alldata[mat]['HOMO'].append(HOMO)
            alldata[mat]['LUMO'].append(LUMO)
            alldata[mat]['iHOMO'].append(iHOMO)
            alldata[mat]['iLUMO'].append(iLUMO)
            alldata[mat]['gap'].append(eg)
            alldata[mat]['fermi'].append((HOMO+LUMO)/2)
            alldata[mat]['eigs']['up'].append(tempeigsup)
            if alldata[mat]['Nspin'] == 2:
                alldata[mat]['eigs']['dn'].append(tempeigsdn)

    ydata = np.array(allgaps[mat])
    #ydata = ydata - ydata[-1]
    line1, = ax.plot(range(len(basissets)),ydata,marker='o',zorder=20,clip_on=False,color=colors[imat])
 #   line1, = ax.plot(allbasis[mat],ydata,marker='o',zorder=20,clip_on=False)
    lines.append(line1)

    #ax.set_xlim([0,6])
    ax.set_ylim([6,15])
    #ax.axhline(y=0,color='k',linewidth=1.5,zorder=10)
    #ax.axvline(x=0,color='k',linewidth=1,zorder=2)
    #if basissets[f] == ['aug-cc-pVQZ','aug-cc-pVQZ']:
    #    ax.set_title(mattitles[mat] + ', ' + titles[f] + ', E$_g$ = {:4.2f} eV'.format(eg) + '\n' + 'Ti and H basis: ' + basissets[f][0],y=1.03)
    #else:
    #    ax.set_title(mattitles[mat] + ', ' + titles[f] + ', E$_g$ = {:4.2f} eV'.format(eg) + '\n' + 'Ti basis: ' + basissets[f][0] + ', H basis: ' + basissets[f][1],y=1.03)
    ax.set_xlabel('Basis set')
    ax.set_ylabel('E$_{gap}$ (eV)')
    ax.set_xticks(range(len(basissetlabels)))
    ax.set_xticklabels(basissetlabels,rotation=90)
    #ax.set_yticklabels([])
    #ax.set_yticks([])
    ax.legend(tuple(lines),tuple(matlabels),loc='upper right',ncol=2,fontsize=11,framealpha=0.95 )

    [ax.spines[x].set_linewidth(3) for x in ax.spines]
    [ax.spines[x].set_zorder(10) for x in ax.spines]
    ax.tick_params(direction = 'out',width = 3,length=5)
    ax.grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig('hydride-gaps.png')

writejson('eigdata.json',alldata)


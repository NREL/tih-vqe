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


mat = 'tih'
#mat = 'lih'
#mat = 'nah'
#mat = 'kh'
files = ['ccsd-sto3g','ccsd-321g','ccsd-631g','ccsd-ccd','ccsd-augccd','ccsd-ccq','ccsd-augccq']
#files = ['ccsd-sto3g','ccsd-631g','ccsd-aug']
#files = ['ccsd-sto3g','ccsd-631g','ccsd-aug','b3lyp-sto3g','b3lyp-631g','b3lyp-aug']

eigdata = readjson('eigdata.json')
eigdata = eigdata[mat]

atomtypes = {'tih':['Ti','H'],'kh':['K','H'],'nah':['Na','H'],'lih':['Li','H']}
atomtypes = atomtypes[mat]

colormats = {'tih':['m','dodgerblue','xkcd:emerald','orange','r','k','purple'],
             'kh':['m','dodgerblue','orange'],
             'nah':['m','dodgerblue','orange','k','r','purple','b','g','y'],
             'lih':['m','dodgerblue','orange','k','r','purple','b','g','y']}
colors = colormats[mat]
ylims = {'tih':[-10,10],'kh':[-10,10],'nah':[-10,10],'lih':[-10,10]}
ylim = ylims[mat]
ymin = ylim[0]
ymax = ylim[1]

mults = {'tih':4,'nah':1,'kh':1,'lih':1}
mult = mults[mat]
mattitles = {'tih':'TiH','kh':'KH','nah':'NaH','lih':'LiH'}
titles = {'ccsd-sto3g':'CCSD','ccsd-321g':'CCSD','ccsd-631g':'CCSD','ccsd-6311g':'CCSD','ccsd-ccd':'CCSD',
          'ccsd-augccd':'CCSD','ccsd-ccq':'CCSD','ccsd-augccq':'CCSD','ccsd-631gdp':'CCSD',
          'b3lyp-631g':'B3LYP','b3lyp-aug':'B3LYP','b3lyp-sto3g':'B3LYP'}
basissets = {'ccsd-sto3g':['STO-3G','STO-3G'],'ccsd-631g':['6-31G','STO-3G'],'ccsd-321g':['3-21G','STO-3G'],
             'ccsd-6311g':['6-311G','STO-3G'],'ccsd-631gdp':['6-31(d\',p\')','STO-3G'],
             'ccsd-ccd':['cc-pVDZ','STO-3G'],'ccsd-augccd':['aug-cc-pVDZ','STO-3G'],
             'ccsd-ccq':['cc-pVQZ','STO-3G'],'ccsd-augccq':['aug-cc-pVQZ','STO-3G'],
          'b3lyp-631g':['6-31G','STO-3G'],'b3lyp-sto3g':['STO-3G','STO-3G'],'b3lyp-aug':['aug-cc-pVQZ','STO-3G']}
atomorbs = {'Ti':['s','p','d','f'],'H':['s'],'K':['s','p'],'Na':['s','p'],'Li':['s','p']}
legorder = {'tih':['H s','Ti s','Ti p','Ti d','Ti f'],
            'lih':['H s','Li s','Li p'],
            'nah':['H s','Na s','Na p'],
            'kh':['H s','K s','K p']}
legorder = legorder[mat]
colordict = dict(zip(legorder,colors[:len(legorder)]))

asq = 0.003
a = asq**0.5

ev2ryd = 13.605698066
ryd2ha = 2.0

xlims = {'tih':[-10,12],'kh':[-5,15],'nah':[-5,15],'lih':[-5,15]}
xlim = xlims[mat]
xmin = xlim[0]
xmax = xlim[1]
xs = 0.01
xrng = np.arange(xmin,xmax+xs/2,xs)
Nvals = len(xrng)
pdos = {x:{y:{'up':np.zeros(Nvals)} for y in atomorbs[x]} for x in atomtypes}
pdos['Total'] = {'up':np.zeros(Nvals)}
Nspin = 1
spindict = {1:['up'],2:['up','dn']}
for f in files:
    basis = f.split('-')[1]

    fig = plt.figure(figsize=(4,4),dpi=200)
    ax = fig.add_subplot(111)

    c = 0
    lines = []
    labels = []

    HOMO = -1000000
    LUMO = 1000000

    filename = 'experimental-length/'+mat+'-'+str(mult)+'-'+f
    if not os.path.isfile(filename):
        print(filename+' not found')
        continue

    index = [i for i in range(len(eigdata['basis'])) if eigdata['basis'][i] == basis][0]
    if eigdata['error'][index]:
        print('error with data',f)
        continue
    HOMO = eigdata['HOMO'][index]
    LUMO = eigdata['LUMO'][index]
    eg = eigdata['gap'][index]
    ef = eigdata['fermi'][index]
    Nbasis = eigdata['Nbasis'][index]
    Nspin = eigdata['Nspin']

    shift = ef
    shift = HOMO + a*2

    data = readfile(filename)


    addedtolegend = []
    startline = find_key('Atomic contributions to Alpha',data)+1
    for line in range(startline,startline+Nbasis):
        temp = data[line].split()
        eig = temp[3]
        eig = float(eig.split('=')[1])*ev2ryd*ryd2ha - shift
        eigdos = delta_func(xrng, a, eig)
        pdos['Total']['up'] += eigdos

        if eig > xmin and eig < xmax:
            projterms = temp[5:]
            for projterm in projterms:
                projterm = projterm.split('=')
                proj = float(projterm[1])
                atorb = projterm[0].split('-')
                at = re.sub('[0-9]+', '', atorb[0])
                orb = atorb[1]
                atorbdos = eigdos*proj
                if orb in pdos[at]:
                    pdos[at][orb]['up'] += atorbdos
                    atorb = at+' '+orb
                else:
                    continue

                legflag = False
                if atorb not in addedtolegend and atorb in legorder:
                    addedtolegend.append(atorb)
                    c += 1
                    legflag = True
                if atorb in legorder:
                    line1, = ax.plot(xrng,atorbdos,color=colordict[atorb],linewidth=1,zorder=5-proj)
                if legflag:
                    lines.append(line1)
                    labels.append(atorb)
                if eig <= 0 and atorb in legorder:
                    ax.fill_between(xrng,np.zeros(Nvals),atorbdos,where=None,color=colordict[atorb],alpha=1,zorder=5-proj)

    if Nspin == 2:
        startline = find_key('Atomic contributions to Beta',data)+1
        pdos['Total']['dn'] = np.zeros(Nvals)
        Nspin = 2
        for line in range(startline,startline+Nbasis):
            temp = data[line].split()
            eig = temp[3]
            eig = float(eig.split('=')[1])*ev2ryd*ryd2ha - shift
            eigdos = delta_func(xrng, a, eig)
            pdos['Total']['dn'] += -1*eigdos

            if eig > xmin and eig < xmax:
                projterms = temp[5:]
                for projterm in projterms:
                    projterm = projterm.split('=')
                    proj = float(projterm[1])
                    atorb = projterm[0].split('-')
                    at = re.sub('[0-9]+', '', atorb[0])
                    orb = atorb[1]
                    if 'dn' not in pdos[at][orb]:
                        pdos[at][orb]['dn'] = np.zeros(Nvals)
                    atorbdos = -1*eigdos*proj
                    if orb in pdos[at]:
                        pdos[at][orb]['dn'] += atorbdos
                        atorb = at+' '+orb
                    else:
                        continue

                    legflag = False
                    if atorb not in colordict and atorb in legorder:
                        colordict[atorb] = colors[c]
                        c += 1
                        legflag = True
                    if atorb in legorder:
                        line1, = ax.plot(xrng,atorbdos,color=colordict[atorb],linewidth=1,zorder=5-proj)
                    if legflag:
                        lines.append(line1)
                        labels.append(atorb)
                    if eig <= 0 and atorb in legorder:
                        ax.fill_between(xrng,np.zeros(Nvals),atorbdos,where=None,color=colordict[atorb],alpha=1,zorder=5-proj)

    templeg = list(zip(labels,lines))
    templeg.sort(key=lambda l: legorder.index(l[0]))
    labels,lines = zip(*templeg)

    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.axhline(y=0,color='k',linewidth=1.5,zorder=10)
    ax.axvline(x=0,color='k',linewidth=1,zorder=2)
    ax.set_title(mattitles[mat] + ', ' + titles[f] + ', E$_g$ = {:4.2f} eV'.format(eg) + '\n' + atomtypes[0]+' basis: ' + basissets[f][0],y=1.03)
    ax.set_xlabel('E - E$_{HOMO}$ (eV)')
    ax.set_ylabel('Density of states')
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.legend(tuple(lines),tuple(labels),loc='lower left',ncol=2,fontsize=11,framealpha=0.95 )

    [ax.spines[x].set_linewidth(3) for x in ax.spines]
    [ax.spines[x].set_zorder(30) for x in ax.spines]
    ax.tick_params(direction = 'out',width = 3,length=5)
    ax.grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(filename+'.png')


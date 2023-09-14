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
#mat = 'nah'

eigdata = readjson('eigdata.json')

mats = ['lih','nah','kh','tih']
matlabels = ['LiH','NaH','KH','TiH']
mults = {'lih':[1],'nah':[1],'kh':[1],'tih':[2,4]}
mults = {'lih':[1],'nah':[1],'kh':[1],'tih':[4]}
basissets = ['sto3g','321g','631g','ccd','augccd','ccq','augccq']
basissetlabels = ['STO-3G','3-21G','6-31G','cc-pVDZ','aug-cc-pVDZ','cc-pVQZ','aug-cc-pVQZ']

indexmats = {'tih':{'atorb':{'H s':0,'Ti s':1,'Ti p':2,'Ti d':3},'order':['H s','Ti s','Ti p','Ti d']},
           'lih':{'atorb':{'H s':0,'Li s':1,'Li p':2},'order':['H s','Li s','Li p']},
           'nah':{'atorb':{'H s':0,'Na s':1,'Na p':2},'order':['H s','Na s','Na p']},
           'kh':{'atorb':{'H s':0,'K s':1,'K p':2},'order':['H s','K s','K p']}}

ev2ryd = 13.605698066
ryd2ha = 2.0
asq = 0.003
a = asq**0.5
thr = 0.1

colors = ['m','cornflowerblue','limegreen','sandybrown']
for imat,mat in enumerate(mats):
    lines = []
    fig = plt.figure(figsize=(6,6),dpi=100)
    ax = fig.add_subplot(211)

    indexes = indexmats[mat]['atorb']
    indexorder = indexmats[mat]['order']
    frontierfracs = np.zeros([4,len(basissets)])
    Nfrontierfracs = np.zeros([4,len(basissets)])
    for m in mults[mat]:
        for ib,b in enumerate(basissets):
            data = readfile('experimental-length/'+mat+'-'+str(m)+'-ccsd-'+b)
            if eigdata[mat]['error'][ib]:
                frontierfracs[:,ib] = np.nan
                continue

            Nbasis = eigdata[mat]['Nbasis'][ib]
            Nspin = eigdata[mat]['Nspin']

            LUMO = eigdata[mat]['LUMO'][ib]

            ieig = 0
            isp = 0
            startline = find_key('Atomic contributions to Alpha',data)+1
            for line in range(startline,startline+Nbasis):
                temp = data[line].split()
                eig = temp[3]
                eig = float(eig.split('=')[1])*ev2ryd*ryd2ha
                if [isp,ieig] == eigdata[mat]['iLUMO'][ib] or abs(LUMO-eig) <= thr:
                    projterms = temp[5:]

                    for projterm in projterms:
                        projterm = projterm.split('=')
                        proj = float(projterm[1])
                        atorb = projterm[0].split('-')
                        at = re.sub('[0-9]+', '', atorb[0])
                        orb = atorb[1]
                        atorb = at+' '+orb
                        index = indexes[atorb]
                        frontierfracs[index,ib] += proj
                        Nfrontierfracs[index,ib] += 1
                    ieig += 1
                else:
                    ieig += 1


            if Nspin == 2:
                isp = 1
                ieig = 0
                startline = find_key('Atomic contributions to Beta',data)+1
                for line in range(startline,startline+Nbasis):
                    temp = data[line].split()
                    eig = temp[3]
                    eig = float(eig.split('=')[1])*ev2ryd*ryd2ha
                    if [isp,ieig] == eigdata[mat]['iLUMO'][ib] or abs(LUMO-eig) <= thr:
                        temp = data[line].split()
                        projterms = temp[5:]

                        for projterm in projterms:
                            projterm = projterm.split('=')
                            proj = float(projterm[1])
                            atorb = projterm[0].split('-')
                            at = re.sub('[0-9]+', '', atorb[0])
                            orb = atorb[1]
                            atorb = at+' '+orb
                            index = indexes[atorb]
                            frontierfracs[index,ib] += proj
                            Nfrontierfracs[index,ib] += 1
                        ieig += 1
                    else:
                        ieig += 1

    Nfrontierfracs[Nfrontierfracs == 0] = 1
    for i in range(len(indexorder)):
        io = indexorder[i]
        iio = indexes[io]
        line1, = ax.plot(range(len(basissets)),frontierfracs[iio,:]/Nfrontierfracs[iio,:],marker='o',
                         color=colors[i],zorder=50,clip_on=False)
        lines.append(line1)
    ax.set_ylim([0,1])
    ax.set_xlim([-0.5,len(basissets)-0.5])
    ax.set_ylabel('LUMO\nprojection')
    ax.set_title(matlabels[imat]+ ' CCSD',y=1.08)
    ncol = 3
    if mat == 'tih':
        ncol = 2
    ax.legend(tuple(lines),tuple(indexorder),fontsize=12,ncol=ncol)
    ax.set_xticklabels([])
    ax.set_xticks(range(len(basissets)))

    [ax.spines[x].set_linewidth(3) for x in ax.spines]
    [ax.spines[x].set_zorder(30) for x in ax.spines]
    ax.tick_params(direction = 'out',width = 3,length=5)
    ax.grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax.set_axisbelow(True)


    ax = fig.add_subplot(212)

    frontierfracs = np.zeros([4,len(basissets)])
    Nfrontierfracs = np.zeros([4,len(basissets)])
    for m in mults[mat]:
        for ib,b in enumerate(basissets):
            data = readfile('experimental-length/'+mat+'-'+str(m)+'-ccsd-'+b)
            if eigdata[mat]['error'][ib]:
                frontierfracs[:,ib] = np.nan
                continue

            Nbasis = eigdata[mat]['Nbasis'][ib]
            Nspin = eigdata[mat]['Nspin']

            HOMO = eigdata[mat]['HOMO'][ib]

            ieig = 0
            isp = 0
            startline = find_key('Atomic contributions to Alpha',data)+1
            for line in range(startline,startline+Nbasis):
                temp = data[line].split()
                eig = temp[3]
                eig = float(eig.split('=')[1])*ev2ryd*ryd2ha
                if [isp,ieig] == eigdata[mat]['iHOMO'][ib] or abs(HOMO-eig) <= thr:
                    temp = data[line].split()
                    projterms = temp[5:]
                    for projterm in projterms:
                        projterm = projterm.split('=')
                        proj = float(projterm[1])
                        atorb = projterm[0].split('-')
                        at = re.sub('[0-9]+', '', atorb[0])
                        orb = atorb[1]
                        atorb = at+' '+orb
                        index = indexes[atorb]
                        frontierfracs[index,ib] += proj
                        Nfrontierfracs[index,ib] += 1
                    ieig += 1
                else:
                    ieig += 1

            if Nspin == 2:
                isp = 1
                ieig = 0
                startline = find_key('Atomic contributions to Beta',data)+1
                for line in range(startline,startline+Nbasis):
                    temp = data[line].split()
                    eig = temp[3]
                    eig = float(eig.split('=')[1])*ev2ryd*ryd2ha
                    if [isp,ieig] == eigdata[mat]['iHOMO'][ib] or abs(HOMO-eig) <= thr:
                        temp = data[line].split()
                        projterms = temp[5:]
                        for projterm in projterms:
                            projterm = projterm.split('=')
                            proj = float(projterm[1])
                            atorb = projterm[0].split('-')
                            at = re.sub('[0-9]+', '', atorb[0])
                            orb = atorb[1]
                            atorb = at+' '+orb
                            index = indexes[atorb]
                            frontierfracs[index,ib] += proj
                            Nfrontierfracs[index,ib] += 1
                        ieig += 1
                    else:
                        ieig += 1

    Nfrontierfracs[Nfrontierfracs == 0] = 1
    for i in range(len(indexorder)):
        io = indexorder[i]
        iio = indexes[io]
        line1, = ax.plot(range(len(basissets)),frontierfracs[iio,:]/Nfrontierfracs[iio,:],marker='o',
                         color=colors[i],zorder=50,clip_on=False)
        lines.append(line1)

    ax.set_xlim([-0.5,len(basissets)-0.5])
    ax.set_ylim([0,1])
    ax.set_xticks(range(len(basissets)))
    ax.set_xticklabels(basissetlabels,rotation=90)

    #ax.axhline(y=0,color='k',linewidth=1.5,zorder=10)
    #ax.axvline(x=0,color='k',linewidth=1,zorder=2)
    #if basissets[f] == ['aug-cc-pVQZ','aug-cc-pVQZ']:
    #    ax.set_title(mattitles[mat] + ', ' + titles[f] + ', E$_g$ = {:4.2f} eV'.format(eg) + '\n' + 'Ti and H basis: ' + basissets[f][0],y=1.03)
    #else:
    #    ax.set_title(mattitles[mat] + ', ' + titles[f] + ', E$_g$ = {:4.2f} eV'.format(eg) + '\n' + 'Ti basis: ' + basissets[f][0] + ', H basis: ' + basissets[f][1],y=1.03)
    ax.set_xlabel(matlabels[imat][:-1]+' basis set')
    ax.set_ylabel('HOMO\nprojection')
    #ax.set_yticklabels([])
    #ax.set_yticks([])
    #ax.legend(tuple(lines),tuple(labels),loc='lower left',ncol=2,fontsize=11,framealpha=0.95 )

    [ax.spines[x].set_linewidth(3) for x in ax.spines]
    [ax.spines[x].set_zorder(30) for x in ax.spines]
    ax.tick_params(direction = 'out',width = 3,length=5)
    ax.grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
    ax.set_axisbelow(True)

    plt.tight_layout()
    # plt.savefig('experimental-length/'+mat+'-'+ f+'.png')

    print(mat)
    for x in eigdata[mat]['gap']:
        print(x)


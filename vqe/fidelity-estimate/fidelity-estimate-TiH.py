#!/usr/bin/env python


import os,sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from matplotlib import cm,ticker
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pylab as pl


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


matstoplot = ['H2','TiHsto3g1','TiHsto3g2','TiHsto3g3']
ucctypes = ['UCCSD','UCCS']
for ucc in ucctypes:
    if ucc == 'UCCSD':
        mats = {'H2':{'label':'H$_2$ STO-3G','Nq':2,'Ns':11,'Nd':4,'ticks':[0.99,0.1,0.01,0.001,0.0001,0.00001],'spam':['1e-1','1e-2','1e-3','1e-4']},
                'TiHsto3g1':{'label':'TiH STO-3G','Nq':6,'Ns':482,'Nd':952,'ticks':[0.1,0.01,0.001,0.0001,0.00001,0.000001],'spam':['1e-2','1e-3','1e-4','1e-5']},
                'TiHsto3g2':{'label':'TiH STO-3G','Nq':10,'Ns':2131,'Nd':5692,'ticks':[0.01,0.001,0.0001,0.00001,0.000001,0.0000001],'spam':['1e-3','1e-4','1e-5','1e-6']},
                'TiHsto3g3':{'label':'TiH STO-3G','Nq':14,'Ns':4935,'Nd':16372,'ticks':[0.01,0.001,0.0001,0.00001,0.000001,0.0000001],'spam':['1e-3','1e-4','1e-5','1e-6']}}
    elif ucc == 'UCCS':
        mats = {'H2':{'label':'H$_2$ STO-3G','Nq':2,'Ns':6,'Nd':0.0001,'ticks':[0.99,0.1,0.01,0.001,0.0001,0.00001],'spam':['1e-1','1e-2','1e-3','1e-4']},
                'TiHsto3g1':{'label':'TiH STO-3G','Nq':6,'Ns':50,'Nd':48,'ticks':[0.1,0.01,0.001,0.0001,0.00001,0.000001],'spam':['1e-2','1e-3','1e-4','1e-5']},
                'TiHsto3g2':{'label':'TiH STO-3G','Nq':10,'Ns':114,'Nd':168,'ticks':[0.01,0.001,0.0001,0.00001,0.000001,0.0000001],'spam':['1e-3','1e-4','1e-5','1e-6']},
                'TiHsto3g3':{'label':'TiH STO-3G','Nq':14,'Ns':178,'Nd':352,'ticks':[0.01,0.001,0.0001,0.00001,0.000001,0.0000001],'spam':['1e-3','1e-4','1e-5','1e-6']}}

    dpi = 300

    fspam = ['1e-3']
    for ispam in fspam:
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
        fig = plt.figure(figsize=(6,6),dpi=dpi)
        gs1 = gridspec.GridSpec(1,1)
        axmain = fig.add_subplot(gs1[0],frameon=False)
        Nx = 2
        Ny = 2
        gs2 = gridspec.GridSpec(Nx,Ny,height_ratios=[1 for x in range(Ny)],width_ratios=[1 for x in range(Nx)])

        ax = []  #list of ax
        for i in range(0,Nx*Ny):
            ax.append(fig.add_subplot(gs2[i]))
        gs2.tight_layout(fig,rect=[0.07,0.1,0.86,0.88])
        gs2.update(hspace=0.47)
        gs2.update(wspace=0.16)
        gs1.tight_layout(fig,rect=[0,0.093,0.93,0.885])

        axmain.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
        axmain.set_xticks([])
        axmain.set_yticks([])
        axmain.set_xticklabels([])
        axmain.set_yticklabels([])
        axmain.grid(False)
        axmain.set_zorder(1)
        #axmain.yaxis.labelpad = 5
        #axmain.set_visible(False)
        axmain.set_title('VQE '+ucc+' fidelity error rate dependence, e$_{q}$ = '+ispam,
                             y=1.15,x=0.56,fontsize=16)

        alpha=0.4
        NL = 120
        coolwarm = cm.get_cmap('rainbow',NL)
        coolwarm = truncate_colormap(coolwarm,0.15,0.85)  #set alpha to 0.4

        lvls = np.arange(-10,0,0.02)
        lvls = np.linspace(-50,0,NL)
        lvls = np.linspace(-2,0,NL)
        cticks = np.arange(-50,1,5)
        cticks = np.arange(-2,1,0.5)

        contour = axmain.contourf([0,0],[0,0],[[0,0],[0,0]], lvls,cmap=coolwarm,alpha=alpha)
        divider = make_axes_locatable(axmain)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        #cticks = np.arange(-10,1,2)
        ctickformat = None
        label = 'log(F)'
        cbar = fig.colorbar(contour, cax=cax, label= label,ticks=cticks,format=ctickformat)
        cbar.outline.set_linewidth(0)
        for i,m in enumerate(matstoplot):
            #data
            factor = 1.1
            fsingle = [mats[m]['ticks'][0]]
            fsingle = [0.05]
            j = 0
            while fsingle[-1] > 1e-7:
                fsingle.append(fsingle[0]/factor**(j+1))
                j += 1
            fsingle = np.array(fsingle)
            fst = np.array(mats[m]['ticks'])
            fst = np.array([0.05,0.01,0.001,0.0001,0.00001,0.000001])
            fdouble = [mats[m]['ticks'][0]]
            fdouble = [0.05]
            j = 0
            while fdouble[-1] > 1e-7:
                fdouble.append(fdouble[0]/factor**(j+1))
                j += 1
            fdouble = np.array(fdouble)
            fdt = np.array(mats[m]['ticks'])
            fdt = np.array([0.05,0.01,0.001,0.0001,0.00001,0.000001])
            fspamlabels = [ispam]                #just keep this constant throughout
            fspam = np.array([float(x) for x in fspamlabels])
            fx,fy = np.meshgrid(fsingle,fdouble)

            matlabel = mats[m]['label']
            Nq = mats[m]['Nq']
            Ns = mats[m]['Ns']
            Nd = mats[m]['Nd']

            #actually plot data
            xmin = fst[-1]
            xmax = fst[0]
            xs = 0.005
            ymin = fdt[-1]
            ymax = fdt[0]
            ys = 0.005

            for k,xspam in enumerate(fspam):
                F = (1.0-fx)**Ns*(1.0-fy)**Nd*(1.0-xspam)**Nq
                F = np.log(F)

                ax[i].set_zorder(2)
                contour = ax[i].contourf(fx,fy,F, lvls,cmap=coolwarm,alpha=alpha)

                benchmarks = [0.95,0.80,0.5,0.1]
                blabels = ['95%','80%','50%','10%']
                for bb,b in enumerate(benchmarks):
                    y = 1-(b/(1-xspam)**Nq/(1-fsingle)**Ns)**(1.0/Nd)
                    spot = [fsingle[x] for x in range(len(y)) if y[x] > fdt[-1]]
                    if spot != []:
                        spot = spot[0]
                    else:
                        continue
                    spot = (np.log(spot)-np.log(fst[-1]))/(np.log(fst[0])-np.log(fst[-1]))
                    if i != 0 or bb <= 1:
                        ax[i].plot(fsingle,y,linewidth=1,color='k',linestyle='-')
                        ax[i].annotate(blabels[bb],rotation=90,xy=(0.001,0.001),xytext=(spot*0.93,0.14),
                                   ha='center',va='center',color='k',fontsize=10,xycoords='axes fraction')

                # divider = make_axes_locatable(ax[i])
                # cax = divider.append_axes('right', size='5%', pad=0.1)
                # cticks = np.arange(-6,1,1)
                # ctickformat = None
                # label = 'label'
                # cbar = fig.colorbar(contour, cax=cax, label= label,ticks=cticks,format=ctickformat)
                # cbar.outline.set_linewidth(0)
                lowc = matplotlib.colors.rgb2hex(coolwarm(0))
                ax[i].patch.set_facecolor('gainsboro')   #is effectively np.nan color
                #ax[i].patch.set_facecolor(lowc)   #is effectively np.nan color
                #ax[i].patch.set_alpha(0.85)

                ax[i].set_xscale('log')
                ax[i].set_yscale('log')
                ax[i].minorticks_off()

                ax[i].axhline(y=0, color='k', linestyle='-',zorder=5)
                ax[i].axvline(x=0, color='k', linestyle='-',zorder=5)
                ax[i].set(xlim=[xmin,xmax],ylim=[ymin,ymax])
                ax[i].set_xticks(np.arange(xmin,xmax+xs,xs))
                ax[i].set_xticks(fst)
                xticks = [round(x,3) for x in np.arange(xmin,xmax+xs,xs)]
                xticks = fst
                ax[i].set_yticks(fdt)
                yticks = [round(x,3) for x in np.arange(ymin,ymax+ys,ys)]
                yticks = fdt
                ax[i].set_yticklabels(yticks)
                if i%Ny == 0:
                    ax[i].set_ylabel('e$_{g_2}$')
                    ax[i].set_yticklabels(yticks)
                else:
                    ax[i].set_yticklabels([])
                if i >= Nx*(Ny-1):
                    ax[i].set_xlabel('e$_{g_1}$')
                    ax[i].set_xticklabels(xticks,rotation=90)
                else:
                    ax[i].set_xticklabels([])

                #ax[i].annotate('e$_{SPAM}$='+fspamlabels[k],ha='center',va='center',fontsize=12,xycoords='axes fraction',
                #               xy=(0.0001,0.0001),xytext=(0.7,0.94))
                ax[i].set_title('{:}, N$_q$={:}\nN$_s$={:}, N$_d$={:}'.format(mats[m]['label'],Nq,Ns,Nd),
                                y=1.03,linespacing=0.95)

        for i in range(0,Nx*Ny):
            [ax[i].spines[x].set_linewidth(3) for x in ax[i].spines]
        #    [ax[i].spines[x].set_zorder(3) for x in ax[i].spines]
            ax[i].tick_params(direction = 'out',width = 3,length=5)
            ax[i].tick_params(top='off',right='off')
            #ax[i].grid(color='gainsboro',linestyle='-',zorder = 1,axis='both')
            ax[i].set_axisbelow(True)

        #plt.tight_layout()
        plt.savefig('overall-fidelity-estimation-'+ucc+'-espam'+ispam+'-tih.png')
        plt.show()


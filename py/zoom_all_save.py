from guacho_utils import *
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import struct
from palettable.cubehelix import red_16
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator

#import seaborn as sns
#sns.set(style="white", color_codes=True,font_scale=1.5)
#sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8, 'xtick.direction': u'out','ytick.direction': u'out','font.family': [u'sans-serif']})

plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['{Palatino}']})

mpl.rcParams.update({'font.size': 20})

rc('text', usetex=True)

mpl.rc('xtick',labelsize=16)
mpl.rc('ytick',labelsize=16)
mpl.rc('text',usetex=True)

rho_HI = 1.66e-24

rmin=100.
rmax=1e6

nout = 80
name=['B1.0','B1.1','B1.5','B5.1','B5.5']

plt.ion()
f, axes = plt.subplots(5,3,sharex='col',sharey='row',figsize=(6,10))
f.subplots_adjust(hspace=0,wspace=0)

for i in range(5):
    ax1 = axes[i][0]
    ax2 = axes[i][1]
    ax3 = axes[i][2]
    j = name[i]

    rr=np.linspace(-13.4,-4.7,150)
    pp=np.linspace(-4.35,4.35,150)
    ex=[-13.4,-4.7,-4.35,4.35] 

    X,Y = np.meshgrid(rr[::],pp[::])

    fileout = open('./DAT/HD20-'+j+'dat',"r")
    npzfile = np.load(fileout)
    print npzfile.files
    rho = npzfile['rho']
    y = npzfile['y']
    by = npzfile['by']
    bx = npzfile['bx']
    b = npzfile['b']
    t = npzfile['t']
    fileout.close()

    rho = rho[24:174,0:150]
    y   =   y[24:174,0:150]
    by  =  by[24:174,0:150]
    bx  =  bx[24:174,0:150]
    b   =   b[24:174,0:150]
    t   =   t[24:174,0:150]

    #campo
    G=by 
    L=bx 
    
    #RHO#
    ax1.set_ylabel('$y[R_{\star}]$')
    if(i==5):
        ax1.set_xlabel('$x[R_{\star}]$')
    #if(i==0):
	#ax1.set_title(j)

    cmap = plt.get_cmap('viridis')
    vmin = 1.0e-20/rho_HI
    vmax = 1.0e-18/rho_HI
    im = ax1.imshow(rho/rho_HI,norm=LogNorm(),cmap=cmap,origin='lower',vmin=vmin,vmax=vmax,extent=ex)
    ax1.contour(y,levels=[.85,.95,.999],colors='black',linewidths=0.5,extent=ex)
    if i==4:
        pos = ax1.get_position()
        axcb1 = f.add_axes([pos.x0,pos.y0-0.06,pos.width,0.02])
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
        cb1 = mpl.colorbar.ColorbarBase(axcb1,cmap=cmap,norm=norm,extend='max',extendfrac=None,orientation='horizontal')
        cb1.set_label(r'$\mathrm{n}_{H}~[cm^{-3}]$')
     
    #TEMP#
    #if(i==0):
#	ax2.set_title(j)

    #if(i==5):
    #    ax2.set_xlabel('$x[R_{\star}]$')

    cmap = 'gist_heat'
    vmin = 1.0e4
    vmax = 2.0e6
    ax2.imshow(t,norm=LogNorm(), origin='lower',cmap=cmap,vmin=vmin,vmax=vmax,extent=ex)
    if i==4:
        pos = ax2.get_position()
        axcb2 = f.add_axes([pos.x0,pos.y0-0.06,pos.width,0.02])
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
        cb2 = mpl.colorbar.ColorbarBase(axcb2,cmap=cmap,norm=norm,extend='max',extendfrac=None,orientation='horizontal')
        cb2.set_label(r'$T~[K] $')

    #BFIELD#
    #ax3.set_xlabel('$x[R_{\star}]$')
   # if(i==0):
#	ax3.set_title(j)

    #if(i==5):
    #    ax3.set_xlabel('$x[R_{\star}]$')

    cmap = 'Spectral'
    vmin = 5.0e-5
    vmax = 5.0e0
    ax3.imshow(b,origin='lower',norm=LogNorm(),cmap=cmap,vmin=vmin,vmax=vmax,extent=ex)

    ax3.set_facecolor('black')

    if i==4:
        pos = ax3.get_position()
        axcb3 = f.add_axes([pos.x0,pos.y0-0.06,pos.width,0.02])
        norm = mpl.colors.LogNorm(vmin=vmin,vmax=vmax)
        cb3 = mpl.colorbar.ColorbarBase(axcb3,cmap=cmap,norm=norm,extend='max',extendfrac=None,orientation='horizontal')
        cb3.set_label(r'$B~[G]$')


    ax3.streamplot(rr,pp,L,G,density=[1,1],color='k',arrowsize=0.5,linewidth=0.5,minlength=0.1)
    ax3.set_xlim(rr.min(),rr.max())
    ax3.set_ylim(pp.min(),pp.max())

plt.show()
plt.ioff()
#plt.savefig("zoom_all.pdf")

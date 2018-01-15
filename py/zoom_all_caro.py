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

rc('text', usetex=True)

mpl.rc('xtick',labelsize=10)
mpl.rc('ytick',labelsize=10)
mpl.rc('text',usetex=True)

rmin=100.
rmax=1e6

fig = plt.figure(figsize=(7, 7))
fig.clf()
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(3,3),
                 axes_pad=0.1,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="each",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )


nout = 80
folder = '/datos_diable/carito/Guacho-master/E'
name=['HD20-B1.1','HD20-B1.5','HD20-B5.5']
#f, ((ax11,ax12,ax13,ax14),(ax21,ax22,ax23,ax24),(ax31,ax32,ax33,ax34)) = plt.subplots(3,4,sharex='col',sharey='row',figsize=(10,10))
#f.subplots_adjust(hspace=0,wspace=0)

#axes = [[ax11,ax12,ax13,ax14],[ax21,ax22,ax23,ax24],[ax31,ax32,ax33,ax34]]

for ax, i in zip(grid.axes_column,name):
    ax1 = ax[0]
    ax2 = ax[1]
    ax3 = ax[2]
    j = i
    path=folder+j+'/BIN/'
    

    rhosc = get_scalings(nout=0,path=path,verbose=False)[2]
    nxtot, nytot, nztot = get_boxsize(nout=0,path=path,verbose=False)
    print nxtot,nytot,nztot
    rho   = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=True)
    vx    = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=True)
    vy    = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=True)
    vz    = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=True)
    pgas  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
    bx    = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)
    by    = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=True)
    bz    = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=True)
    rho_n = readbin3d_all(nout=nout,neq=8,path=path,verbose=False,mhd=True)*rhosc
    passives = readbin3d_all(nout=nout,neq=9,path=path,verbose=False,mhd=True)*rhosc
    phi   = readbin3d_all(nout=nout,neq=0,path=path,base='ph-',verbose=False)    

    v=np.sqrt(vx*vx+vy*vy+vz*vz)
    b=np.sqrt(bx*bx+by*by+bz*bz)
    t=pgas*1./(8.3145e7*(2*rho-rho_n))  #para H_RATE
    y=1.-rho_n/rho
    phi=phi#+1e-30
    #psi=(phi*rho_n/rhosc)*2.1789e-11
    psi=(phi*rho_n)*2.1789e-11
    
    rr=np.linspace(-0.15,-0.08,200)
    pp=np.linspace(-0.035,0.035,200)
    ex=[-0.15,-0.08,-0.035,0.035] 

    xi=0
    xf=200
    yi=15
    yf=215
    
    X,Y = np.meshgrid(rr[::],pp[::])

    #campo
    G=by[nztot/2,yi:yf,xi:xf] 
    L=bx[nztot/2,yi:yf,xi:xf] 
    
    #RHO#
    #if(i==0):
    ax1.set_ylabel('y[AU]')
    ax1.set_title(j)

    cmap = red_16.mpl_colormap
    vmin = np.min(psi)#1.0e-20
    vmax = np.max(psi)#1.0e-18
    im = ax1.imshow(psi[nztot/2,yi:yf,xi:xf],cmap=cmap,origin='lower',extent=ex,norm=LogNorm())
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cb1=ax1.cax.colorbar(im,norm=norm)
    ax1.cax.toggle_label(True)
    cb1.set_label_text(r'$\rho~[gcm^{-3}]$')
     
    #TEMP#
   # if(i==0):
    ax2.set_ylabel('y[AU]')

    cmap = 'gist_heat'
    vmin = 1.e4
    vmax = 1.e5
    im2=ax2.imshow(t[nztot/2,yi:yf,xi:xf],norm=LogNorm(), origin='lower',cmap=cmap,vmin=vmin,vmax=vmax,extent=ex)
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cb2=ax2.cax.colorbar(im2,norm=norm)
    ax2.cax.toggle_label(True)
    cb2.set_label_text(r'$T~[K]$')

    #BFIELD#
    ax3.set_xlabel('x[AU]')
    #if(i==0):
    ax3.set_ylabel('y[AU]')

    cmap = 'Spectral'
    vmin =10.
    vmax =240
    im3=ax3.imshow(v[nztot/2,yi:yf,xi:xf]*1e-5,origin='lower',cmap=cmap,extent=ex,vmax=vmax,vmin=vmin,norm=LogNorm())
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cb3=ax3.cax.colorbar(im3,norm=norm)
    ax3.cax.toggle_label(True)
    cb3.set_label_text(r'$Psi~$')

    #ax3.streamplot(X,Y,L,G,density=[2,2],color='k')
    #ax3.set_xlim(rr.min(),rr.max())
    #ax3.set_ylim(pp.min(),pp.max())

plt.savefig("./PDF/zoom_all_two.pdf")
plt.show()

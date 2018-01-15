from guacho_utils import *
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import struct
from palettable.cubehelix import red_16

#import seaborn as sns
#sns.set(style="white", color_codes=True,font_scale=1.5)
#sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8, 'xtick.direction': u'out','ytick.direction': u'out','font.family': [u'sans-serif']})

plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['{Palatino}']})

rc('text', usetex=True)

mpl.rc('xtick',labelsize=10)
mpl.rc('ytick',labelsize=10)
mpl.rc('text',usetex=True)

plt.ion()
rmin=100.
rmax=1e6

nout = 45
run=['carito/Guacho-master/EHD20-B1.5','carito/Guacho-master/EHD20-B5.5','carito/Guacho-master/EHD20-B5.1','matias/Guacho_Working/Guacho-master/EHD20-B1.1']
name=['HD20-B1.5','HD20-B5.5','HD20-B5.1','HD20-B1.1']

#for k in range(len(run)):
for k in range(0,4):
    i = run[k]
    j = name[k]
    path ='/datos_diable/'+i+'/BIN/'
    
    rhosc = get_scalings(nout=0, path=path, verbose=True )[2]
    tempsc = get_scalings(nout=0, path=path, verbose=True )[2]
    nxtot, nytot, nztot = get_boxsize(nout=0, path=path, verbose=True )
    rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=True)
    vx = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=True)
    vy = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=True)
    vz = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=True)
    pgas = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
    bx = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)
    by = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=True)
    bz = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=True)
    rho_n = readbin3d_all(nout=nout,neq=8,path=path,verbose=False,mhd=True)*rhosc
    phi  = readbin3d_all(nout=nout,neq=0,path=path,base='ph-',verbose=True)
    print nxtot,nytot,nztot
    
    v=np.sqrt(vx*vx+vy*vy+vz*vz)
    b=np.sqrt(bx*bx+by*by+bz*bz)
    t=pgas*1./(8.3145e7*(2*rho-rho_n))  #para H_RATE
    y=1.-rho_n/rho
    phi=phi + 1e-30
    psi=(phi*rho_n/rhosc)*2.1789e-11
    #pt=pgas+b*b/(8*np.pi) +rho*v*v
    #va=b/np.sqrt(4*np.pi*rho)
    #vs=np.sqrt(1.01*pgas/rho)
    #ca=v/va

    if (k >= 2):
       rr=np.linspace(0.08,0.15,182)
       pp=np.linspace(-0.035,0.035,182)
       ex=[0.08,0.15,-0.035,0.035]
       #vale para la salida 45
       zi=218
       zf=400
       yi=8
       yf=190

    else:

       rr=np.linspace(0.08,0.15,200)
       pp=np.linspace(-0.035,0.035,200)
       ex=[0.08,0.15,-0.035,0.035] 
       #vale para la salida 45
       zi=260
       zf=460
       yi=15
       yf=215
    

    Z,Y = np.meshgrid(rr[::],pp[::])
    
    
    #velocidad
    U=vy[zi:zf,yi:yf,nxtot/2].transpose()
    W=vz[zi:zf,yi:yf,nxtot/2].transpose()
    #campo
    G=by[zi:zf,yi:yf,nxtot/2].transpose() 
    L=bz[zi:zf,yi:yf,nxtot/2].transpose() 
    
    Z = Z#.transpose()
    Y = Y#.transpose()
    
    plt.figure(1)
    plt.clf()
    plt.title(j)
    plt.xlabel('z[AU]')
    plt.ylabel('y[AU]')
    plt.imshow(v[zi:zf,yi:yf,nxtot/2].transpose()/1e5, origin='lower', cmap='rainbow',extent=ex,vmin=0,vmax=400)
    cb= plt.colorbar(extend='max')
    cb.set_label(r'$V~[km~s^{-1}]$')
    #plt.quiver(Z[::3],Y[::3],W[::3],U[::3],pivot='mid',units='x')
    #plt.savefig('./PDF/'+j+'-V.pdf',dpi=300,transparent=True,bbox_inches='tight')
    
    #
    
    plt.figure(2)
    plt.clf()
    plt.title(j)
    plt.xlabel('z[AU]')
    plt.ylabel('y[AU]')

    plt.imshow(b[zi:zf,yi:yf,nxtot/2].transpose(), origin='lower',norm=LogNorm(), cmap='Spectral',extent=ex,vmin=5e-5, vmax=5e0)
    cb=plt.colorbar(extend='max')
    cb.set_label(r'$B~[G]$')
    plt.streamplot(Z,Y,L,G,density=[2,2],color='k')
    plt.xlim(rr.min(),rr.max())
    plt.ylim(pp.min(),pp.max())
    #plt.savefig('./PDF/'+j+'-B.pdf',dpi=300,transparent=True,bbox_inches='tight')
    ###
    ####
    ###
    plt.figure(3)
    plt.clf()
    plt.title(j)
    plt.xlabel('z[AU]')
    plt.ylabel('y[AU]')

    plt.imshow(t[zi:zf,yi:yf,nxtot/2].transpose(),norm=LogNorm(), origin='lower',cmap='gist_heat',extent=ex,vmax=1.5e6)
    cb=plt.colorbar(extend='max')
    cb.set_label(r'$T~[K]$')
    #plt.savefig('./PDF/'+j+'-T.pdf',dpi=300,transparent=True,bbox_inches='tight')
    ###


    #cmap = sns.cubehelix_palette(light=1, as_cmap=True)
    cmap = red_16.mpl_colormap
    
    plt.figure(4)
    plt.clf()
    plt.title(j)
    plt.xlabel('z[AU]')
    plt.ylabel('y[AU]')

    plt.imshow(rho[zi:zf,yi:yf,nxtot/2].transpose(),norm=LogNorm(),cmap=cmap,origin='lower',extent=ex,vmin=1.e-20,vmax=1.e-18)
    cb=plt.colorbar(extend='max')
    plt.contour(y[zi:zf,yi:yf,nxtot/2].transpose() ,levels=[.85, .95, .999],colors='black',lw=1,extent=ex)
    cb.set_label(r'$\rho~[gcm^{-3}]$')
    #plt.savefig('./PDF/'+j+'-rho.pdf',dpi=300,transparent=True,bbox_inches='tight')

    plt.figure(5)
    plt.clf()
    plt.title(j)
    plt.xlabel('z[AU]')
    plt.ylabel('y[AU]')


    plt.imshow(rho_n[zi:zf,yi:yf,nxtot/2].transpose(),norm=LogNorm(), origin='lower',cmap='gist_heat',extent=ex)#,vmax=1.5e6)
    cb=plt.colorbar(extend='max')
    cb.set_label(r'$rho_n~[K]$')
    plt.savefig('./PDF/'+j+'-rho_n.pdf',dpi=300,transparent=True,bbox_inches='tight')
    ###





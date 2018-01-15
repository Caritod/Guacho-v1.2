#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import * 
from matplotlib.colors import LogNorm

plt.ion()
rmin=100.
rmax=1e6

#name=['EHD20-B1.1','EHD20-B1.5','EHD20-B5.1','EHD20-B5.5']
name=['SW']
nout =29

lx=13.4  #para hd20
ly=lx/2.   #
lz=lx   #

for i in range(0,1):
    
    #path='/datos_diable/carito/Guacho-master/'+name[i]+'/BIN/'
    path='../EXO/'+name[i]+'/BIN/'
    rhosc = get_scalings(nout=0, path=path, verbose=True )[2]
    nxtot, nytot, nztot = get_boxsize(nout=0, path=path, verbose=True )
    print nxtot
    rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=True)
    vx = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=True)
    vy = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=True)
    vz = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=True)
    pgas = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
    bx = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)
    by = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=True)
    bz = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=True)
    rho_n = readbin3d_all(nout=nout,neq=8,path=path,verbose=False,mhd=True)*rhosc
    #phi  = readbin3d_all(nout=nout,neq=0,path=path,base='ph-',verbose=True)
        
    v=np.sqrt(vx*vx+vy*vy+vz*vz)
    b=np.sqrt(bx*bx+by*by+bz*bz)
    t=pgas*1./(8.3145e7*(2.*rho-rho_n))
    #t=pgas*0.5*1.66e-24/(1.38e-16*rho)
    va=b/np.sqrt(4*np.pi*rho)
    ca=v/va
    vs=np.sqrt(2.*1.38e-16*t/1.66e-24)
    cs=v/vs
    beta=8*np.pi*pgas/b**2
    #phi=phi#+1e-30
    #psi=(phi*rho_n/rhosc)*2.1789e-11
    #y=1.-rho_n/rho
    
    XX = np.linspace(-lx,lx,nxtot)
    ZZ = np.linspace(-lz,lz,nztot)
    YY = np.linspace(-ly,ly,nytot)
    ex=[-lx,lx,-lz,lz] #en radio solar
    ey=[-lx,lx,-ly,ly] #en radio solar
    
    plt.figure(1)
    plt.clf()
    plt.title(r'Total Density')
    plt.xlabel(r'$x[R_\odot]$')
    plt.ylabel(r'$z[R_\odot]$')
    plt.imshow(rho_n[nztot/2,::,::], norm=LogNorm(),origin='lower', cmap='Blues')#,vmin= 1e-20, vmax=1e-17,interpolation='none')
    plt.colorbar()
    #plt.savefig('./PDF/'+name[i]+'-rhoXZ.pdf',dpi=300,transparent=True,bbox_inches='tight')
    print('****')
    
    plt.figure(2)
    plt.clf()
    nskip=4
    #X,Z = np.meshgrid(XX[::nskip],ZZ[::nskip])
    X,Y = np.meshgrid(XX[::nskip],YY[::nskip])
    #U=vx[::nskip,nytot/2,::nskip]
    U=bx[nztot/2,::nskip,::nskip]
    #W=vz[::nskip,nytot/2,::nskip]
    W=by[nztot/2,::nskip,::nskip]
        
    plt.title(r'B')
    plt.xlabel(r'$x[R_\odot]$')
    plt.ylabel(r'$z[R_\odot]$')
    plt.imshow(b[nztot/2,::,::], origin='lower', cmap='rainbow',extent=ey,norm=LogNorm())#,vmax=500 , vmin=1)
    plt.colorbar()
    plt.contour(cs[nztot/2,::,::],colors='white', vmin=1, levels=[1],extent=ey,lw=3,ls='--' )
    plt.contour(ca[nztot/2,::,::],colors='black', vmin=1, levels=[1],extent=ey,lw=3,ls='--' )
    #plt.quiver(X,Y,U,W,pivot='mid',units='x')#,scale=1e12)
    plt.streamplot(X,Y,U,W,density=[3,3],color='k')
   # plt.savefig('./PDF/'+name[i]+'-BXY-62.pdf',dpi=300,transparent=True,bbox_inches='tight')
    
    print('****')
    plt.figure(3)
    plt.clf()
    
    plt.title(r'Temperature')
    plt.xlabel(r'$x[R_\odot]$')
    plt.ylabel(r'$y[R_\odot]$')
    plt.imshow(t[nztot/2,::,::],origin='lower',cmap='gist_heat',norm=LogNorm())
    plt.colorbar()
    #plt.savefig('./PDF/'+name[i]+'TXZ.pdf',dpi=300,transparent=True,bbox_inches='tight')
    print('****')

    plt.figure(4)
    plt.clf()
    plt.title(r'Velocity')
    plt.xlabel(r'$x[R_\odot]$')
    plt.ylabel(r'$y[R_\odot]$')
    
    nskip=4
    X,Z = np.meshgrid(XX[::nskip],ZZ[::nskip])
    W=vz[::nskip,nytot/2,::nskip]
    U=vx[::nskip,nytot/2,::nskip]
    plt.imshow(v[::,nytot/2,::]/1e5,origin='lower',cmap='rainbow',extent=ex)
    plt.colorbar()

    plt.quiver(X,Z,U,W,pivot='mid',units='x')#,scale=1e12)
    #plt.savefig('./PDF/'+name[i]+'TXZ.pdf',dpi=300,transparent=True,bbox_inches='tight')
    print('****')

#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import * 
from matplotlib.colors import LogNorm

plt.ion()
rmin=100.
rmax=1e6

path = '/datos_diable/carito/Guacho-master/EHD20-B5.5/BIN/'
nout = 13


rhosc = get_scalings(nout=0, path=path, verbose=True )[2]
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
v=np.sqrt(vx*vx+vy*vy+vz*vz)
b=np.sqrt(bx*bx+by*by+bz*bz)
t=pgas*1./(8.3145e7*(2.*rho-rho_n))
#t=pgas*0.5*1.66e-24/(1.38e-16*rho)
va=b/np.sqrt(4*np.pi*rho)
ca=v/va
vs=np.sqrt(2.*1.38e-16*1.5e6/1.66e-24)
cs=v/vs
beta=8*np.pi*pgas/b**2
phi=phi#+1e-30
psi=(phi*rho_n/rhosc)*2.1789e-11
y=1.-rho_n/rho
lx=11.472  #para nx#tot=256
ly=5.7359#6.7218   #
lz=11.472#13.444   #
#lx=7.15 #para nxtot=100
#ly=lx#4.12
#lz=lx
XX = np.linspace(-lx,lx,nxtot)
ZZ = np.linspace(-lz,lz,nztot)
YY = np.linspace(-ly,ly,nytot)
#ex=[-5.3774,5.3774,-5.3774,5.3774] #en radio solar
ex=[-lx,lx,-lz,lz] #en radio solar
ey=[-lx,lx,-ly,ly] #en radio solar

plt.figure(1)
plt.clf()
plt.title(r'Total Density')
plt.xlabel(r'$x[R_\odot]$')
plt.ylabel(r'$z[R_\odot]$')
plt.imshow(rho[::,nytot/2,::], norm=LogNorm(),origin='lower', cmap='Blues',extent=ex)
plt.colorbar()
#plt.savefig('./PDF/20b1-D2-Dxz.pdf',dpi=300,transparent=True,bbox_inches='tight')
print('****')

plt.figure(2)
plt.clf()
nskip=4
X,Z = np.meshgrid(XX[::nskip],ZZ[::nskip])
U=vx[::nskip,nytot/2,::nskip]
W=vz[::nskip,nytot/2,::nskip]

plt.title(r'Velocity')
plt.xlabel(r'$x[R_\odot]$')
plt.ylabel(r'$z[R_\odot]$')
plt.imshow(v[::,nytot/2,::]/1e5, origin='lower', cmap='rainbow',extent=ex)#,vmax=500 )
plt.colorbar()
plt.contour(cs[::,nytot/2,::],colors='white', vmin=1, levels=[1],extent=ex,lw=3,ls='--' )
plt.quiver(X,Z,U,W,pivot='mid',units='x')
#plt.savefig('./PDF/20b1-D2-Vxz.pdf',dpi=300,transparent=True,bbox_inches='tight')
#plt.xlabel(r'$x~[\mathrm{UA}]$')
#plt.ylabel(r'$z~[\mathrm{UA}]$')

print('****')
plt.figure(3)
plt.clf()

plt.title(r'Temperature')
plt.xlabel(r'$x[R_\odot]$')
plt.ylabel(r'$y[R_\odot]$')
plt.imshow(t[::,nytot/2,::],origin='lower',cmap='gist_heat',extent=ex)#max=2e6)
plt.colorbar()
#plt.savefig('./PDF/20b1-Txy.pdf',dpi=300,transparent=True,bbox_inches='tight')
print('****')


plt.figure(4)
plt.clf()
nskip=4
K,Y = np.meshgrid(XX[::nskip],YY[::nskip])
M=by[nztot/2,::nskip,::nskip]
N=bx[nztot/2,::nskip,::nskip]
plt.title(r'Magnetic Field')
plt.xlabel(r'$x[R_\odot]$')
plt.ylabel(r'$y[R_\odot]$')

plt.imshow(b[nztot/2,::,::],origin='lower',cmap='Blues',norm=LogNorm(),extent=ey)#,vmin=1e-3 )
plt.colorbar()
#plt.savefig('./PDF/20b1-D2-Bxy.pdf',dpi=300,transparent=True,bbox_inches='tight')
plt.streamplot(K,Y,N,M,density=[2,2],color='k')
print('****')

#x_ax, y_ax, z_ax = get_axis(nout=0, path=path,verbose=False)
plt.figure(5)
plt.clf()
nskip=4
K,Y = np.meshgrid(XX[::nskip],YY[::nskip])
H=vx[nztot/2,::nskip,::nskip]
J=vy[nztot/2,::nskip,::nskip]

plt.title(r'Velocity')
plt.xlabel(r'$x[R_\odot]$')
plt.ylabel(r'$y[R_\odot]$')
plt.imshow(v[nztot/2,::,::]/1e5, origin='lower', cmap='rainbow',extent=ey)#,vmax=1800, vmin=150)
plt.colorbar()
plt.contour(ca[nztot/2,::,::],colors='white', vmin=1, levels=[1],extent=ey,lw=3)
plt.contour(cs[nztot/2,::,::],colors='white', vmin=1, levels=[1],extent=ey,lw=3,ls='--' )
#plt.streamplot(K,Y,N,M,density=[2,2],color='k')
#plt.savefig('./PDF/20b1-D2-Vyz.pdf',dpi=300,transparent=True,bbox_inches='tight')
#plt.quiver(K,Y,H,J,pivot='mid',units='width')

print('****')
#yph=1.-rho_n/rho
plt.figure(6)
plt.clf()

plt.title(r'psi')
plt.xlabel(r'$x[R_\odot]$')
plt.ylabel(r'$y[R_\odot]$')

plt.imshow(psi[::,nytot/2,::], origin='lower',norm=LogNorm(),cmap='Spectral',extent=ex )#vmin=1e-18 , vmax= 1e-12 
plt.colorbar()

plt.figure(7)
plt.imshow(y[::,nytot/2,::], origin='lower', cmap='OrRd',extent=ex)
plt.colorbar()
#plt.savefig('./PDF/20b1-Pxy.pdf',dpi=300,transparent=True,bbox_inches='tight')
print('****')



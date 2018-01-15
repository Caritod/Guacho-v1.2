#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib as mpl

#plt.ion()
fig = plt.figure(figsize=(7, 7))
fig.clf()
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(2,3),
                 axes_pad=0.5,
                 share_all=False,
                 #aspect=True,
                 cbar_location="right",
                 cbar_mode="each",
                 cbar_size="7%",
                 cbar_pad=0,
                 )


rmin=100.
rmax=1e6
Rsun=6.955e10
lx= 1.
ly= 3.
ex= [-lx,lx,0,ly]

path = '../prom2d/BIN/'
nout = 2

nxtot, nytot, nztot = get_boxsize(nout=0, path=path, verbose=True )
rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=True)
pgas = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
t=pgas*1.66e-24*0.5/(1.38e-16*rho)
bz = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=True)
by = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=True)
bx = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)
b=np.sqrt(bx**2+by**2+bz**2)
vz = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=True)
vy = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=True)
vx = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=True)
v=np.sqrt(vx**2+vy**2+vz**2)

#for i in zip(grid.axes_column):
#    ax1 = ax[0]
#    ax2 = ax[1]
#    ax3 = ax[2]
#    ax4 = ax[3]
#    ax5 = ax[4]
#    ax6 = ax[6]

grid[0].set_ylabel('y[Rs]')
grid[0].set_title(r'$\rho$ [g/cm3]')
grid[0].set_xlabel('x[Rs]')
vmin=np.min(rho)
vmax=np.max(rho)
im0=grid[0].imshow(rho[nztot/2,::,::], origin='lower', cmap='Blues',extent=ex,vmin=vmin,vmax=vmax )
norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
cb1=grid[0].cax.colorbar(im0,norm=norm)
grid[0].cax.toggle_label(True)
#cb1.set_label_text(r'$\beta $')

grid[1].set_ylabel('y[Rs]')
grid[1].set_title('T [K]')
grid[1].set_xlabel('x[Rs]')
im1=grid[1].imshow(t[nztot/2,::,::], origin='lower', cmap='gist_heat',extent=ex )
cb2=grid[1].cax.colorbar(im1)
#grid[1].cax.toggle_label(True)
#cb2.set_label_text(r'$\beta $')

XX = np.linspace(-1.,1.,nxtot)
#ZZ = np.linspace(-lz,lz,nztot)
YY = np.linspace(0,3.,nytot)
nskip=6
X,Y = np.meshgrid(XX[::nskip],YY[::nskip])
U=bx[nztot/2,::nskip,::nskip]
W=by[nztot/2,::nskip,::nskip]

grid[2].set_ylabel('y[Rs]')
grid[2].set_title('B [G]')
grid[2].set_xlabel('x[Rs]')
im2=grid[2].imshow(b[nztot/2,::,::], origin='lower', cmap='Blues',extent=ex)
cb3=grid[2].cax.colorbar(im2)
grid[2].cax.toggle_label(True)
#cb3.set_label_text(r'$\beta $')
grid[2].streamplot(X,Y,U,W,density=[1,1],color='k')

T=vx[nztot/2,::nskip,::nskip]
P=vy[nztot/2,::nskip,::nskip]
grid[3].set_ylabel('y[Rs]')
grid[3].set_xlabel('x[Rs]')
grid[3].set_title('V [km/s]')
im3=grid[3].imshow(v[nztot/2,::,::]/1e5, origin='lower', cmap='rainbow' ,extent=ex)
cb4=grid[3].cax.colorbar(im3)
grid[3].cax.toggle_label(True)
#cb4.set_label_text(r'$\beta $')
grid[3].quiver(X,Y,T,P,color='k',pivot='mid',units='xy',width=.01,scale_units='inches',scale=1e9)

beta=pgas*8*np.pi/b**2
grid[4].set_ylabel('y[Rs]')
grid[4].set_xlabel('x[Rs]')
grid[4].set_title(r'$\beta$')
im4=grid[4].imshow(beta[nztot/2,::,::], origin='lower', cmap='rainbow' ,extent=ex,vmax=1)
cb5=grid[4].cax.colorbar(im4)
grid[4].cax.toggle_label(True)
#cb5.set_label_text(r'$\beta $')


grid[5].set_ylabel('y[Rs]')
grid[5].set_title('B [G]')
grid[5].set_xlabel('x[Rs]')
im2=grid[5].imshow(t[::,nytot/2,::], origin='lower', cmap='gist_heat',extent=ex)
cb3=grid[5].cax.colorbar(im2)
grid[5].cax.toggle_label(True)
#cb3.set_label_text(r'$\beta $')
#grid[5].streamplot(X,Y,U,W,density=[1,1],color='k')
plt.show()
#plt.ioff()

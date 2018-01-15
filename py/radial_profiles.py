import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm

amh=1.66e-24
rs=8.3460e+10
Rjup=7.1492E9
rp=1.38*Rjup
au=1.496e13
dx=0.15*au/460   #(xphys=0.15AU)
x=(np.arange(0,230,dtype=float)+0.5)*dx/rp
a=0.047*au/rs

dx1=0.15*au/400   #(xphys=0.15AU)
x1=(np.arange(0,200,dtype=float)+0.5)*dx/rp

fileout0=open('DAT/EHD20-B1.0_1D-80.txt.npz','r')
fileout05=open('DAT/EHD20-B1.05_1D-80.txt.npz','r')
fileout1=open('DAT/EHD20-B1.1_1D-80.txt.npz','r')
fileout2=open('DAT/EHD20-B1.5_1D-80.txt.npz','r')
fileout3=open('DAT/EHD20-B5.5_1D-80.txt.npz','r')
fileout4=open('DAT/EHD20-B5.1_1D-80.txt.npz','r')
data00=np.load(fileout0)
data05=np.load(fileout05)
data11=np.load(fileout1)
data15=np.load(fileout2)
data55=np.load(fileout3)
data51=np.load(fileout4)

dd=np.array([data00['d'],data05['d'], data11['d'], data15['d'], data55['d'],data51['d']])
dd=dd/amh
tt=np.array([data00['t'],data05['t'], data11['t'], data15['t'], data55['t'],data51['t']])
vv=np.array([data00['v'],data05['v'], data11['v'], data15['v'], data55['v'],data51['v']])
bb=np.array([data00['b'],data05['b'], data11['b'], data15['b'], data55['b'],data51['b']])
ppg=np.array([data00['pg'],data05['pg'], data11['pg'], data15['pg'], data55['pg'],data51['pg']])

pd=0.5*dd*vv**2*amh
pm=bb**2/(8.*np.pi)

#####label = ['B1.0','B1.05','B1.1','B1.5','B5.5','B5.1']
#####
#####fig = plt.figure(figsize=(6, 4))
#####grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
#####                 nrows_ncols=(4,1),
#####                 axes_pad=0.1,
#####                 aspect=True,
#####                 share_all=False
#####                 )
#####
#####indx=np.where((x>14.7) & (x<24.7))[0]
#####
#####
######for i in range(4):
######    grid[0].plot(x[indx]-11.7,np.log10(dd[i,indx]),label=label[i])
######    grid[1].plot(x[indx]-11.7,np.log10(bb[i,indx]),label=label[i])
######    grid[2].plot(x[indx]-11.7,np.log10(vv[i,indx]),label=label[i])
######    grid[3].plot(x[indx]-11.7,np.log10(tt[i,indx]),label=label[i])
#####
#####
#####
#####
######grid[2].plot(x[indx]-11.7,np.log10(pm[1,indx]),label='Pmag')
######grid[2].plot(x[indx]-11.7,np.log10(ppg[1,indx]),label='Pgas')
######grid[2].plot(x[indx]-11.7,np.log10(pd[1,indx])-20,label='Pdyn')
#####
#####grid[0].set_ylim(4,7)
#####grid[0].set_ylabel(r'$\rho$')
#####grid[2].set_ylabel(r'V')
#####grid[1].set_ylim(-5,0)
#####grid[1].set_ylabel(r'B')
#####grid[2].set_ylim(6,8)
#####grid[2].set_xlabel(r'R $[R_p]$')
#####
#####grid[2].set_xlim(3,8)
#####
#####
######grid[1].set_yticks([-2, 0, 2])
######grid[0].legend()
#####grid[1].legend()
######grid[2].legend()
#####
######plt.savefig('./PDF/radial_profiles.pdf',transparent=True,bbox_inches='tight',dpi=300)




fig,axs = plt.subplots(1,1,figsize=(15,5))

label = ['B1.1','B1.5','B5.1']

axs.set_title(label[2])
axs.plot(x-42.0,ppg[0],label='pg')
axs.plot(x-42.0, pd[0],label='pd')
axs.plot(x-42.0, pm[0],label='pm')
axs.plot(x-42.0,ppg[2],'r',label='pg')
axs.plot(x-42.0, pd[2],'r',label='pd')
axs.plot(x-42.0, pm[2],'r',label='pm')

#axs.plot(x1-36.24,ppg[5],label='pg')
#axs.plot(x1-36.24, pd[5],label='pd')
#axs.plot(x1-36.24, pm[5],label='pm')
#q=ind[i]
#ax.set_title(k)
#ax.set_xlabel(r'$x~[\mathrm{R_*}]$')
#ax.set_ylabel(r'$Pressure$')

axs.legend()
axs.set_yscale('log')
plt.show()

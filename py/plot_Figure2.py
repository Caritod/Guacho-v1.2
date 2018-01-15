from guacho_utils import *
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm
from matplotlib import rc
import matplotlib.pyplot as plt
#import yt
plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['{Palatino}']})

rc('text', usetex=True)

mpl.rc('xtick',labelsize=10)
mpl.rc('ytick',labelsize=10)
mpl.rc('text',usetex=True)

firstrun= True
plt.ion()

path='/datos_diable/matias/Guacho_Working/Guacho-master/A1-5GP-GC/BIN/'
nout=45
lx=11.472
ly=5.7359
extent=(-lx,lx,-lx,lx)

circ1=plt.Circle((0.,0.),radius=1.,color='k',fill=False,ls='dashed')
#circ2=plt.Circle((0.,0.),radius=0.035,color='k',fill=False,ls='dashed')
  
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
    

v=np.sqrt(vx*vx+vy*vy+vz*vz)
b=np.sqrt(bx*bx+by*by+bz*bz)
t=pgas*1./(8.3145e7*(2*rho-rho_n))      
yhp=1.-rho_n/rho
#Phi=Phi+1e-30
#Psi= (Phi*rho_n)*2.1789e-11

#------------------------------------------------------ 
fig= plt.figure(1)
fig.clf()

im = plt.imshow(rho[:,nytot/2,:], norm=LogNorm(), origin='lower', cmap='Blues',
           extent=extent )

cont= plt.contour(yhp[:,nytot/2,:], levels=[.9, .99, .999], colors='k',extent=extent )

plt.xlabel(r'$x~[\mathrm{R_*}]$')
plt.ylabel(r'$z~[\mathrm{R_*}]$')
plt.title(r'Density')
#plt.text(-0.13, 0.13,'A1 @ t='+str(nout*0.025).format('f')+' days')
 
ax= fig.add_subplot(1,1,1)
ax.add_patch(circ1)
#ax.add_patch(circ2)

cb= plt.colorbar(im,label='', pad=0.01, extend='max')
cb.set_label(r'$\rho_\mathrm{H}~[\mathrm{g cm^{-3}}]$')

plt.savefig('./PDF/A1b-f1.pdf',dpi=300,transparent=True,bbox_inches='tight')

#-------------------------------------------------------
fig= plt.figure(2)
fig.clf()

im = plt.imshow(t[:,nytot/2,:], origin='lower', cmap='gist_heat',
     extent=extent, vmin= 1e4, vmax = 2e6, norm= LogNorm())

plt.xlabel(r'$x~[\mathrm{R_*}]$')
plt.ylabel(r'$z~[\mathrm{R_*}]$')

plt.title(r'Temperature')

cb= plt.colorbar(im,label='', pad=0.01, extend='max')
cb.set_label(r'$T~[\mathrm{K}]$')
#plt.text(-0.13, 0.13,'P4c @ t='+str(nout*0.025).format('f')+' days')

plt.savefig('./PDF/A1b-f2.pdf',dpi=300,transparent=True,bbox_inches='tight')
#-------------------------------------------------------
fig= plt.figure(3)
fig.clf()

im = plt.imshow(b[:,nytot/2,:], origin='lower', cmap='Spectral', #bds_highcontrast',
     extent=extent, norm= LogNorm())#, vmin=1e-18 , vmax= 1e-12 )
plt.xlabel(r'$x~[\mathrm{R_*}]$')
plt.ylabel(r'$z~[\mathrm{R_*}]$')
plt.title(r'Magnetic field')
#plt.text(-0.13, 0.13,'P4c @ t='+str(nout*0.025).format('f')+' days')

cb= plt.colorbar(im,label='', pad=0.01, extend='max')
#cb.set_label(r'$\psi~[\mathrm{erg~cm^{-3}~s^{-1}}]$')
cb.set_label(r'$B~[G]$')

plt.savefig('./PDF/A1b-f3.pdf',dpi=300,transparent=True,bbox_inches='tight')

#--------------------------------------------------------------------------
fig= plt.figure(4)
fig.clf()

im = plt.imshow(v[:,nytot/2,:]/1e5,origin='lower', cmap='rainbow', extent=extent )

plt.xlabel(r'$x~[\mathrm{R_*}]$')
plt.ylabel(r'$z~[\mathrm{R_*}]$')
plt.title(r'Velocity')
#plt.text(-0.13, 0.13,'P4c @ t='+str(nout*0.025).format('f')+' days')
 
cb= plt.colorbar(im,label='', pad=0.01, extend='max')
cb.set_label(r'$v~[\mathrm{km s^{-1}}]$')

plt.savefig('./PDF/A1b-f4.pdf',dpi=300,transparent=True,bbox_inches='tight')



'''
    data = dict(density = (denT, "g/cm**3"), temperature=(temp, "K"),ionization_fraction= (yhp, "K") )
    bbox = np.array([[-.15, .15], [-.0375, .0375] ,[-.15, .15]])
    ds = yt.load_uniform_grid(data, denT.shape, length_unit="au", bbox=bbox, nprocs=1)

  slc = yt.SlicePlot(ds, "y", ["density"], width=((.3,'au'),(0.3,'au')) )

  slc.set_cmap("density", "Blues_r")
  slc.set_cmap("temperature", "gist_heat")
  slc.set_xlabel('x[AU]')
  slc.set_ylabel('z[AU]')
  slc.set_colorbar_label("density","$n_\mathrm{H}~[\mathrm{cm^{-3}}]$")

  slc.set_zlim("density",zmax=1e6,zmin=1e2)
  
  slc.annotate_sphere([0.,0.,0.,],(0.00511, 'au'),circle_args={'ls': 'solid'})
  slc.annotate_sphere([0.,0.,0.,],(0.035, 'au'),circle_args={'ls': 'dashed'})
  slc.annotate_contour('ionization_fraction',ncont=1,clim=(0.5,0.5))
  slc.annotate_contour('ionization_fraction',ncont=1,clim=(0.99,0.99))
  slc.annotate_contour('ionization_fraction',ncont=1,clim=(0.999,0.999),plot_args={"colors": 'black',"ls":'.'})

  slc.save('dens.png')

'''
  

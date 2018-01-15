import fortranfile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from guacho_utils import *
from mpl_toolkits.axes_grid1 import ImageGrid
#import seaborn as sns
#import yt
plt.rc('font', family='serif')
plt.rc('font',**{'family':'serif','serif':['{Palatino}']})

mpl.rc('xtick',labelsize=10)
mpl.rc('ytick',labelsize=10)
mpl.rc('text',usetex=True)
endian='<'      # little endian
kind='d'        # double precision

velinit=-300
velfinal=300
nvmap=250
dvel=(velfinal-velinit)/float(nvmap)

sintheta=np.sin(-3.33*np.pi/180.)
costheta=np.cos(-3.33*np.pi/180.)

rstar= 1.2*6.955e10
rplan= 1.35*7.1492E9 /rstar
Rorb = 0.047*1.496e13/rstar
torb = 3.52*86400.

plt.ion()
#sns.set_style('ticks')
#sns.set_context('paper', font_scale=1.5, rc={"lines.linewidth": 1.5})
#sns.set_style(rc={'font.family': ['Serif'], 'font.sans-serif':['Palatino']})

fig = plt.figure(figsize=(14.625, 4.5))
fig.clf()
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,3),
                 axes_pad=0.1,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.15,
                 )

nout = 80
folder= '/datos_diable/carito/Guacho-master/E'
name=['HD20-B1.1','HD20-B1.5','HD20-B5.5']

for ax, k in zip(grid,name):

    path=folder+k+'/BIN/'
   
    if (k =='HD20-B5.1'):
      nxmap = 400
      nymap = 200
    else: 
      nxmap = 460
      nymap = 230

    dx = 0.150*1.49e13/float(nxmap) #  dx in cm ( 1AU = 1.49e13 cm)

    rs = 1.2*6.955e10#/1.49e13
    rp = 1.35*7.1492E9 /rs

    time = 0.025*(nout)*86400. 
    xp = Rorb*np.cos(2.*np.pi*time/torb-25.*np.pi/180.)
    zp = Rorb*np.sin(2.*np.pi*time/torb-25.*np.pi/180.)
    yp = xp*sintheta
    
    filein = path+'LA_tau_planeta-'+str(nout).zfill(3)+'.bin'
    f = fortranfile.FortranFile(filein,endian=endian)
    data = f.readReals(kind).reshape(nxmap,nymap,nvmap,order='F').T

    emtaunu = np.exp(-data)
    datap = 1.0-np.sum(emtaunu,0)/250.
    
    ax.set_title(k)
    ax.set_xlabel(r'$z~[\mathrm{R_*}]$')
    ax.set_ylabel(r'$y~[\mathrm{R_*}]$')
   
    axs= fig.add_subplot(ax)
    circ1=plt.Circle((0.,0.),radius=1.,color='w',fill=False)
    circ2=plt.Circle((zp,yp),radius=rplan,color='k',fill=True)

    extent=[-2.9,2.9,-2.9,2.9]

    im= ax.imshow(datap[(nymap/2-51):(nymap/2+49),(nxmap/2-51):(nxmap/2+49)],
               origin='lower',vmin=0.,vmax=.6,extent=extent,cmap='pink')

    axs.add_patch(circ1)
    axs.add_patch(circ2)

    cb=ax.cax.colorbar(im)
    ax.cax.toggle_label(True)
    cb.set_label_text(r'Ly~$\alpha$ absorption fraction $(1-\mathrm{e}^{-\tau})$')


    fileout = './PDF/tau_all.pdf'# % my_name
    plt.savefig(fileout,dpi=300,transparent=True,bbox_inches='tight')

plt.show()


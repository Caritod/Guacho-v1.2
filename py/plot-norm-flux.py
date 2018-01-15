import fortranfile
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


B1_0  = np.genfromtxt('./DAT/EHD20-B1.0-abs-080.dat')
#B1_05 = np.genfromtxt('./DAT/EHD20-B1.05-abs-080.dat')
B1_1  = np.genfromtxt('./DAT/EHD20-B1.1-abs-080.dat')
B1_5  = np.genfromtxt('./DAT/EHD20-B1.5-abs-080.dat')
B5_1  = np.genfromtxt('./DAT/EHD20-B5.1-abs-045.dat')
B5_5  = np.genfromtxt('./DAT/EHD20-B5.5-abs-080.dat')

v_helios=-15.
disc_block=0.985

# Interpolation
vmin =-300
vmax = 300
nvel = 250
x=np.linspace(vmin,vmax,nvel)

#---------plot---------------
plt.ion()
sns.set_style('ticks')
sns.set_context('paper', font_scale=1.5, rc={"lines.linewidth": 1.5})
sns.set_style(rc={'font.family': ['Serif'], 'font.sans-serif':['Palatino']})
plt.figure(figsize=(8, 6))
plt.clf() 

figure_title=''
plt.text(.5, 1.02, figure_title,
         horizontalalignment='center',
         fontsize=17)

plt.xrange=[]
plt.ylim(0.3, 1.0)
plt.xlim(-300,300)
plt.xlabel(r'Velocity [km s$^{-1}$]')
plt.ylabel(r'Normalized Flux')

ls1='--'        #m_p=2
ls2='-'         #m_p=1
c1='#e41a1c'#'#d95f02'   #sns.xkcd_rgb["pale red"]       #f=1x
c2='#4daf4a'#'#1b9e77'    #sns.xkcd_rgb["medium green"]   #f=0.2x
c3='#984ea3'#'#7570b3'   #sns.xkcd_rgb["denim blue"]     #f=5x
c4='#377eb8'#'#e7298a'
['#e41a1c','#377eb8','#4daf4a','#984ea3']
#['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
plt.plot(x-v_helios,B1_0[:,1]*disc_block,c=c2,label='B1.0',ls='-')
#plt.plot(x-v_helios,B1_05[:,1]*disc_block,c=c2,label='B1.05',ls='-')
plt.plot(x-v_helios,B1_1[:,1]*disc_block,c=c1,label='B1.1',ls='-')
plt.plot(x-v_helios,B1_5[:,1]*disc_block,c=c4,label='B1.5',ls='-')
plt.plot(x-v_helios,B5_1[:,1]*disc_block,c=c1,label='B5.1',ls='--')
plt.plot(x-v_helios,B5_5[:,1]*disc_block,c=c4,label='B5.5',ls='--')
plt.grid(True)

plt.axvspan(-40., 40., facecolor='#EACF00', alpha=0.5,label='Geo. abs.')
plt.legend(loc=4,frameon=True,numpoints= 1, prop={'size':12}, ncol=1)

plt.savefig('./PDF/abs-prof.pdf',transparent=True,bbox_inches='tight',dpi=300)


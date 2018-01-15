import fortranfile
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rc('xtick',labelsize=10)
mpl.rc('ytick',labelsize=10)
#mpl.rc("legend", fontsize=40)


data1   = np.genfromtxt('./DAT/in_transit.dat')
data  = np.genfromtxt('./DAT/out_off_transit.dat')
B1_0  = np.genfromtxt('./DAT/EHD20-B1.0-abs-080.dat')
#B1_05 = np.genfromtxt('./DAT/EHD20-B1.05-abs-080.dat')
B1_1  = np.genfromtxt('./DAT/EHD20-B1.1-abs-080.dat')
B1_5  = np.genfromtxt('./DAT/EHD20-B1.5-abs-080.dat')
B5_1  = np.genfromtxt('./DAT/EHD20-B5.1-abs-045.dat')
B5_5  = np.genfromtxt('./DAT/EHD20-B5.5-abs-080.dat')

vel= data[:,0]
EM = data[:,1]
vel1= data1[:,0]
EM1 = data1[:,1]

v_helios=-15.

disc_block=0.985

# Interpolation
vmin =-300
vmax = 300
nvel=250
x=np.linspace(vmin,vmax,nvel)

f = interp1d(vel, EM, kind='linear',fill_value=0)
f1 = interp1d(vel1, EM1, kind='linear',fill_value=0)
#########################################################

plt.ion()
sns.set_style('ticks')
sns.set_context('paper', font_scale=1.5, rc={"lines.linewidth": 1.5})
sns.set_style(rc={'font.family': ['Serif'], 'font.sans-serif':['Palatino']})
#sns.set(font_scale=3)  # crazy big
plt.figure(figsize=(8, 6))
plt.clf() 
#figure_title = r"(a) $v_*= 130 $km s$^{-1}$" #, $\dot{m}_{P}=[1-2]\times 10^{10}$g s$^{-1}$"
#figure_title = r"(b) $v_*= 205 $km s$^{-1}$" #, $\dot{m}_{P}=[1-2]\times 10^{10}$g s$^{-1}$"
#figure_title = r"(c) $v_*= 372 $km s$^{-1}$" #,$\dot{m}_{P}=[1-2]\times 10^{10}$g s$^{-1}$"
figure_title=''
plt.text(.5, 3.75, figure_title,
         horizontalalignment='center',
         fontsize=17)

plt.xrange=[]
plt.ylim(-0.4,3.7)
plt.xlim(-300,300)
plt.xlabel(r'Velocity [km s$^{-1}$]')
plt.ylabel(r'Flux ($\times  10^{14}$ erg cm$^{-2}$ s$^{-1}$ \AA )')

plt.plot(x,f1(x),'-',label=r'In transit (VM03)',lw=2, c='k')
plt.plot(x,f(x),'-',label=r'Out of transit (VM03)',lw=1,c='k' )

ls1='--'        #m_p=2
ls2='-'         #m_p=1
c1='#e41a1c'#'#d95f02'   #sns.xkcd_rgb["pale red"]       #f=1x
c2='#4daf4a'#'#1b9e77'    #sns.xkcd_rgb["medium green"]   #f=0.2x
c3='#984ea3'#'#7570b3'   #sns.xkcd_rgb["denim blue"]     #f=5x
c4='#377eb8'#'#e7298a'


plt.plot(x-v_helios,f(x-v_helios)*B1_0[:,1]*disc_block,c=c2,label='B1.0',ls='-')
#plt.plot(x-v_helios,f(x-v_helios)*B1_05[:,1]*disc_block,c=c2,label='B1.05',ls='-')
plt.plot(x-v_helios,f(x-v_helios)*B1_1[:,1]*disc_block,c=c1,label='B1.1',ls='-')
plt.plot(x-v_helios,f(x-v_helios)*B1_5[:,1]*disc_block,c=c4,label='B1.5',ls='-')
plt.plot(x-v_helios,f(x-v_helios)*B5_1[:,1]*disc_block,c=c1,label='B5.1',ls='--')
plt.plot(x-v_helios,f(x-v_helios)*B5_5[:,1]*disc_block,c=c4,label='B5.5',ls='--')

plt.grid(True)

plt.axvspan(-40., 40., facecolor='#EACF00', alpha=0.5,label='Geo. abs.')
plt.legend(loc='upper left',frameon=True,numpoints = 1, prop={'size':12}, ncol=3)

#plt.savefig('../PNG/NUEVOS/f6-b.png',transparent=True,dpi=300,bbox_inches='tight')
plt.savefig('./PDF/VM-comparison.pdf',transparent=True,bbox_inches='tight',dpi=300)
#plt.savefig('../EPS/NUEVOS/f6-b.eps',transparent=True,bbox_inches='tight',dpi=300)

plt.ioff()

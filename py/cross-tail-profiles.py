from guacho_utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np 
import pylab
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
    
#f1='DAT/EHD20-B1.1.txt'
#f2='DAT/EHD20-B1.5.txt'
#f3='DAT/EHD20-B5.1.txt'
#f4='DAT/EHD20-B5.5.txt'
#---------------------------grafico----------------------------------
data1   = np.genfromtxt('./center-hA1-5GP-GC.txt')
d1b= data1[:,0]
t1b= data1[:,1]
v1b= data1[:,2]

data2   = np.genfromtxt('./center-hA1-DR-GC.txt')
d1= data1[:,0]
t1= data1[:,1]
v1= data1[:,2]

data3   = np.genfromtxt('./center-hA2-5GP-GC.txt')
d2b= data3[:,0]
t2b= data3[:,1]
v2b= data3[:,2]

data4   = np.genfromtxt('./center-hA2-DR-GC.txt')
d2= data4[:,0]
t2= data4[:,1]
v2= data4[:,2]

rs=8.3460e+10
au=1.496e13
dx=0.128/256   #(xphys=0.128AU)
x=np.arange(-44,36)*dx*au/rs

RSW3   = 1.2*rs/au
rp=0.047*au/rs

mk1='ko'
mk2='k<'
lwd=1.3
msz=5

c1='#e41a1c'  #'#d95f02'   #sns.xkcd_rgb["pale red"]       #f=1x
c2='#377eb8'  #'#1b9e77'    #sns.xkcd_rgb["medium green"]   #f=0.2x
c3='#4daf4a'  #'#7570b3'   #sns.xkcd_rgb["denim blue"]     #f=5x
c4='#984ea3'  #'#e7298a'
c5='#ff7f00'  #'#66a61e'

sns.set_style('ticks')
sns.set_context('paper', font_scale=1.5, rc={"lines.linewidth": 1.5})
sns.set_style(rc={'font.family': ['Serif'], 'font.sans-serif':['Palatino']})

plt.ion()
f=plt.figure(1,figsize=(15, 5))
#figure_title = r"Stellar wind parameters"
#plt.text(0.5, 1.12, figure_title,
#         horizontalalignment='center',
#         fontsize=15,
#         transform = ax2.transAxes)

ax1 = f.add_subplot(1, 3, 1)
ax2 = f.add_subplot(1, 3, 2)#,sharex=ax1)
ax3 = f.add_subplot(1, 3, 3)#,sharex=ax1)

ax1.set_title('Density [g cm$^{-3}$]')
ax1.semilogx(d1b,x,c=c1,label='A1b')  
ax1.semilogx(d2b,x,c=c2,label='A2b')  
ax1.set_ylabel(r"r [R$_*$]") 
ax1.legend(loc='upper right',frameon=True,numpoints = 1, prop={'size':10}, ncol=2)

ax2.set_title(r"Temperature [x10$^6$ K]") 
ax2.plot(t1b/1e6,x,c=c1)  
ax2.plot(t2b/1e6,x,c=c2)  
#ax2.axvline(x=rp,ls='-.',c='k')
#ax2.set_ylim(0.2,1.8)
ax2.set_ylabel(r"r [R$_*$]") 

ax3.set_title('Velocity [km s$^{-1}$]')
ax3.plot(v1b/1.e5,x,c=c1)  
ax3.plot(v2b/1.e5,x,c=c2)  
#ax3.axvline(x=rp,ls='-.',c='k')
ax3.set_ylabel(r"r [R$_*$]") 

f.subplots_adjust(hspace=0.3, top=0.85)
#plt.setp(ax1.get_xticklabels(), visible=False)
#plt.setp(ax2.get_xticklabels(), visible=False)
#plt.savefig('../PDF/fig1.png',transparent=True,bbox_inches='tight',dpi=300)
#plt.savefig('./PDF/spl-b-param.pdf',transparent=True,bbox_inches='tight',dpi=300)
#f.set_size_inches(3.,9.)

################################################presiones##########################


pres4   = np.genfromtxt('./center-pA1-DR-GC.txt')
pd4= pres4[:,0]
pm4= pres4[:,1]
pg4= pres4[:,2]
pt4= pd4+pm4+pg4
beta4=pg4/pm4
pres5   = np.genfromtxt('./center-pA2-DR-GC.txt')
pd5= pres5[:,0]
pm5= pres5[:,1]
pg5= pres5[:,2]
pt5= pd5+pm5+pg5
beta5=pg5/pm5
###

r=plt.figure(2,figsize=(10, 5))
ax1 = r.add_subplot(1, 2, 1)
ax2 = r.add_subplot(1, 2, 2,sharey=ax1)

ax1.set_title('A1 model')
ax1.semilogy(x,pd4,c=c1,label='$P_d$')  
ax1.semilogy(x,pm4,c=c2,label='$P_m$')  
ax1.semilogy(x,pg4,c=c3,label='$P_g$')  
ax1.semilogy(x,pt4,c='k',label='$P_t$', lw=lwd)  
#ax1.axvline(x=rp,ls='-.',c='k',label='R$_o$')
ax1.set_xlabel(r"y [R$_*$]") 
#ax1.set_ylim(1e-8,1e2)
ax1.legend(loc='upper right',frameon=True,numpoints = 1, prop={'size':10}, ncol=2)

ax2.set_title(r"A2 model") 
ax2.semilogy(x,pd5,c=c1,label='$P_d$')  
ax2.semilogy(x,pm5,c=c2,label='$P_m$')  
ax2.semilogy(x,pg5,c=c3,label='$P_g$')  
ax2.semilogy(x,pt5,c='k',label='$P_t$',lw=lwd)  
#ax2.axvline(x=rp,ls='-.',c='k',label='R$_o$')
#ax2.set_ylim(1e-8,1e2)
ax2.set_xlabel(r"y [R$_*$]") 

r.subplots_adjust(hspace=0.3, top=0.85)
#plt.savefig('./PDF/zoom/ztail-pres.pdf',transparent=True,bbox_inches='tight',dpi=300)
#

##beta
#pres4   = np.genfromtxt('./wind-presWA1-256-GC.txt')
#pres1   = np.genfromtxt('./spl-presA1-DR-GC.txt')
#pm1= pres1[:,1]
#pg1= pres1[:,2]
#beta1=pg1/pm1
##pres5   = np.genfromtxt('./wind-presWA2-256-GC.txt')
#pres2   = np.genfromtxt('./spl-presA2-DR-GC.txt')
#pm2= pres2[:,1]
#pg2= pres2[:,2]
#beta2=pg2/pm2

plt.figure(3)
plt.clf()
plt.semilogy(x,pm4,c='g',label='$A1$')  
plt.semilogy(x,pm5,c=c2,label='$A2$') 

#f=plt.figure(3,figsize=(5, 5))
#ax1 = f.add_subplot(1, 1, 1)
#
#ax1.set_title(r"Beta")
#ax1.plot(x,beta1,c=c1,label= '$A1$')
#ax1.plot(x,beta2,c=c2,label= '$A2$')
#ax1.plot(x,beta4,c=c3,label= '$A1b$')
#ax1.plot(x,beta5,c=c4,label= '$A2b$')
#ax1.set_xlabel(r"r [R$_*$]")
#ax1.set_ylim(0,1)
#ax1.legend(loc='upper right',frameon=True,numpoints = 1, prop={'size':10}, ncol=2)
#f.subplots_adjust(hspace=0.3, top=0.85)
 

plt.ioff()

from guacho_utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 
import pylab
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def get_WD_analytic():
  
  filename2 = './HD20-B1-WD.dat'
  
  # get number of radial steps in file
  nLines = sum(1 for line in open(filename2))
  nR = nLines - 7
  
  # make arrays
  rGridWD = np.zeros(nR)
  densGridWD = np.zeros(nR)
  VrGridWD = np.zeros(nR)
  VphiGridWD = np.zeros(nR)
  BrGridWD = np.zeros(nR)
  BphiGridWD = np.zeros(nR)
  
  # open file
  WDfile = open(filename2,'r')

  # skip header
  next(WDfile)
  next(WDfile)
  next(WDfile)
  next(WDfile)
  next(WDfile)
  next(WDfile)
  next(WDfile)
  
  # go through and fill arrays
  iR = 0
  for line in WDfile:
    data = line.split()
    
    rGridWD[iR] = data[0]
    densGridWD[iR] = data[1]
    VrGridWD[iR] = data[2]
    VphiGridWD[iR] = data[3]
    BrGridWD[iR] = data[4]
    BphiGridWD[iR] = data[5]
    
    iR += 1
  
  # close file 
  WDfile.close()
  
  return rGridWD,densGridWD,VrGridWD,VphiGridWD,BrGridWD,BphiGridWD

#---------------------------grafico----------------------------------
data1n   = np.genfromtxt('./HD20-Pol1-MHD.txt')
d1n= data1n[:,0]
#t1n= data1n[:,1]
v1n= data1n[:,2]
#b1n= data1n[:,3]

data3a  = np.genfromtxt('./EHD20-B1.1.txt')
d3n= data3a[:,0]
v3n= data3a[:,2]


data5n  = np.genfromtxt('./test20B1-parker.dat')
r5n= data5n[:,0]
d5n= data5n[:,1]
v5n= data5n[:,2]

# get results from anayltic model
rGridWD,densGridWD,VrGridWD,VphiGridWD,BrGridWD,BphiGridWD = get_WD_analytic()
V=np.sqrt(VrGridWD**2+VphiGridWD**2)
B=np.sqrt(BrGridWD**2+BphiGridWD**2)

rs=8.3460e+10
au=1.496e13
dx2=0.05/100   #(xphys=0.128AU)
dx1=0.15/400   #(xphys=0.128AU)
dx3=0.08/200
x1=np.arange(0,200)*dx1*au/rs
x2=np.arange(0,50)*dx2*au/rs
x3=np.arange(0,100)*dx3*au/rs

RSW3   = 1.2*rs/au
rp=0.047*au/rs

plt.figure()

#plt.plot(x2,v2n/1e5,'b')
plt.plot(x3,v1n/1e5,'r')
plt.plot(x1,v3n/1e5,'b')
plt.plot(r5n,v5n,'o')
plt.plot(rGridWD,V,'k')

plt.figure()
plt.semilogy(x3,d1n,'r')
#plt.semilogy(x1,d3n,'g')
plt.semilogy(r5n,d5n,'o')
plt.semilogy(x1,d3n,'b')
plt.semilogy(rGridWD,densGridWD,'k')

#plt.figure()
#plt.plot(x1,b1n,'r')
#plt.plot(x2,b2n,'b')
#plt.plot(rGridWD,B,'k')
plt.show()

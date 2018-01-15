from guacho_utils import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np 
import pylab

firstrun=True
#plt.ion()
runs=['EHD20-B1.1','EHD20-B1.5','EHD20-B5.5','EHD20-B5.1']
nout = 80

for ii in range( 0,3):
    run=runs[ii]
    path = '/datos_diable/carito/Guacho-master/'+run+'/BIN/'

    rhosc = get_scalings(nout=0, path=path, verbose=False )[2]
    tempsc = get_scalings(nout=0, path=path, verbose=False )[2]
    nxtot, nytot, nztot = get_boxsize(nout=0, path=path, verbose=False)
    print 'nxtot=', nxtot
    rho   = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=True)
    vx    = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=True)
    vy    = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=True)
    vz    = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=True)
    pg    = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=True)
    bx  = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)
    by  = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=True)
    bz  = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=True)
    rho_n = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=True)*rhosc

    v=np.sqrt(vx*vx+vy*vy+vz*vz)
    b=np.sqrt(bx*bx+by*by+bz*bz)
    t=pg*1./(8.3145e7*(2*rho-rho_n))  #para H_RATE
#    t=pg*0.5*1.66e-24/(1.38e-16*rho)  #para ADIABATIC

#    va=b/np.sqrt(4*np.pi*rho)
#    vs=np.sqrt(1.1*pgas/rho)
#    ca=v/va
#    cs=v/vs

#    pd=rho*v**2
#    pm=(b**2/8*np.pi)

    d =   rho[nztot/2, nytot/2, 0:nxtot/2]
    t =    t [nztot/2, nytot/2, 0:nxtot/2]
    v =    v [nztot/2, nytot/2, 0:nxtot/2]
    b =    b [nztot/2, nytot/2, 0:nxtot/2]
    rn= rho_n[nztot/2, nytot/2, 0:nxtot/2]

    a=np.array([d,t,v,b,rn])
#    b=np.array([pd,pg])
    np.savetxt('DAT/'+run+'cp-80.txt', np.transpose(a))
#    np.savetxt(run+'.txt', np.transpose(b))


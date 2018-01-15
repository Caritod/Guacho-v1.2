from pyevtk.vtk import VtkFile, VtkRectilinearGrid
import numpy as np
import struct
from guacho_utils import *
import sys 

ifile = int(sys.argv[1])

mhd = True
nskip = 1
nout = 80
folder = '/datos_diable/carito/Guacho-master/E'
names = ['HD20-B1.5','HD20-B5.5','HD20-B5.1','HD20-B1.1']
path = folder + names[ifile] + '/BIN/'
OutputFile = names[ifile]

rhosc = get_scalings(nout=nout, path=path, verbose=True )[2]
nxtot, nytot, nztot = get_boxsize(nout=nout, path=path, verbose=True)
x_extent, y_extent, z_extent = get_extent(nout,path=path,verbose=False)
dx, dy, dz = get_deltas(nout,path=path,verbose=False)

print dx,dy,dz
print x_extent, y_extent, z_extent
print nxtot,nytot,nztot

rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False,mhd=mhd)
vx  = readbin3d_all(nout=nout,neq=1,path=path,verbose=False,mhd=mhd)/1.e5
vy  = readbin3d_all(nout=nout,neq=2,path=path,verbose=False,mhd=mhd)/1.e5
vz  = readbin3d_all(nout=nout,neq=3,path=path,verbose=False,mhd=mhd)/1.e5
pgas  = readbin3d_all(nout=nout,neq=4,path=path,verbose=False,mhd=mhd)
if mhd:
  bx  = readbin3d_all(nout=nout,neq=5,path=path,verbose=False,mhd=mhd)
  by  = readbin3d_all(nout=nout,neq=6,path=path,verbose=False,mhd=mhd)
  bz  = readbin3d_all(nout=nout,neq=7,path=path,verbose=False,mhd=mhd)
rho_n  = readbin3d_all(nout=nout,neq=8,path=path,verbose=False,mhd=mhd)*rhosc


##if ifile == 3:
##  ixmin = 135
##  ixmax = 265
##  iymin = 34
##  iymax = 164
##  izmin = 270
##  izmax = 400
##  nxtot = 130; nytot = 130; nztot = 130
##else:
##  ixmin = 155
##  ixmax = 305
##  iymin = 40
##  iymax = 190
##  izmin = 310
##  izmax = 460
##  nxtot = 150; nytot = 150; nztot = 150

ixmin = 0
ixmax = 460
iymin = 69
iymax = 161
izmin = 0
izmax = 460
nxtot = 460; nytot = 92; nztot = 460

rrho = rho[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
vvx = vx[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
vvy = vy[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
vvz = vz[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
ppgas = pgas[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
if mhd:
  bbx = bx[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
  bby = by[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
  bbz = bz[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]
rrho_n = rho_n[izmin:izmax:nskip,iymin:iymax:nskip,ixmin:ixmax:nskip]

rrho = rrho.astype(dtype="float64")
vvx = vvx.astype(dtype="float64")
vvy = vvy.astype(dtype="float64")
vvz = vvz.astype(dtype="float64")
ppgas = ppgas.astype(dtype="float64")
if mhd:
  bbx = bbx.astype(dtype="float64")
  bby = bby.astype(dtype="float64")
  bbz = bbz.astype(dtype="float64")
rrho_n = rrho_n.astype(dtype="float64")

Vel = np.sqrt(vvx*vvx+vvy*vvy+vvz*vvz)
print 'Vel: %f %f' % (Vel.min(),Vel.max())
print nxtot, nytot, nztot, rhosc
print rho.shape

nx, ny, nz = nxtot/nskip, nytot/nskip, nztot/nskip
lx, ly, lz = 1.0, 0.20, 1.0
dx, dy, dz = lx/nx, ly/ny, lz/nz
ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1)
x = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
z = np.arange(0, lz + 0.1*dz, dz, dtype='float64')
start, end = (0,0,0), (nx, ny, nz)

w = VtkFile(OutputFile, VtkRectilinearGrid)
w.openGrid(start = start, end = end)
w.openPiece( start = start, end = end)

# Point data
#passive = np.zeros(npoints, dtype="float64", order='F')
#w.openData("Point", scalars = "Passive")
#w.addData("Passive", passive)
#w.closeData("Point")

# Cell data
#pressure = np.zeros([nx, ny, nz], dtype="float64", order='F')
if mhd:
  w.openData("Cell", scalars = ("Density","RhoNeutros","GasPressure"), vectors = ("Velocity","Magnetic"))
else:
  w.openData("Cell", scalars = "Density", vectors = ("Velocity"))
w.addData("Density", rrho)
w.addData("Velocity", (vvx,vvy,vvz))
w.addData("GasPressure",ppgas)
if mhd:
  w.addData("Magnetic", (bbx,bby,bbz))
w.addData("RhoNeutros", rrho_n)
w.closeData("Cell")

# Coordinates of cell vertices
w.openElement("Coordinates")
w.addData("x_coordinates", x);
w.addData("y_coordinates", y);
w.addData("z_coordinates", z);
w.closeElement("Coordinates");

w.closePiece()
w.closeGrid()

#w.appendData(data = passive)
w.appendData(data = rrho)
w.appendData(data = (vvx,vvy,vvz))
w.appendData(data = ppgas)
if mhd:
  w.appendData(data = (bbx,bby,bbz))
w.appendData(data = rrho_n)
w.appendData(x).appendData(y).appendData(z)
w.save()
#

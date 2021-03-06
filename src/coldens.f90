!=======================================================================
!> @file coldens.f90
!> @brief Column density projection
!> @author Alejandro Esquivel
!> @date 4/May/2016

! Copyright (c) 2016 Guacho Co-Op
!
! This file is part of Guacho-3D.
!
! Guacho-3D is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================

!> @brief Column density projection
!> @details Utilities to compute a column density map

module coldens_utilities

contains

!> @brief Initializes data
!> @details Initializes data, MPI and other stuff

subroutine init_coldens()

!  Initializes MPI, data arrays, etc
use parameters
use globals, only : u, dx, dy, dz, coords, rank, left, right   &
                     , top, bottom, out, in, rank, comm3d
implicit none
  integer :: nps, err
  integer, dimension(0:ndim-1) :: dims
  logical, dimension(0:ndim-1) :: period

  !initializes MPI
#ifdef MPIP
#ifdef PERIODX
  logical, parameter :: perx=.true.
#else
  logical, parameter :: perx=.false.
#endif
#ifdef PERIODY
  logical, parameter :: pery=.true.
#else
  logical, parameter :: pery=.false.
#endif
#ifdef PERIODZ
  logical, parameter :: perz=.true.
#else
  logical, parameter :: perz=.false.
#endif
  period(0)=perx
  period(1)=pery
  period(2)=perz
  dims(0)  =MPI_NBX
  dims(1)  =MPI_NBY
  dims(2)  =MPI_NBZ
  
  call mpi_init (err)
  call mpi_comm_rank (mpi_comm_world,rank,err)
  call mpi_comm_size (mpi_comm_world,nps,err)
  if (nps.ne.np) then
     print*, 'processor number (',nps,') is not equal to pre-defined number (',np,')'
     call mpi_finalize(err) 
     stop
  endif
#else
  rank=0
  coords(:)=0
#endif
  if(rank.eq.master) then
     print '(a)' ,"*******************************************"
     print '(a)' ,"                        _                 *"
     print '(a)' ,"  __   _   _  __ _  ___| |__   ___    3   *"
     print '(a)' ," / _ `| | | |/ _` |/ __| '_ \ / _ \    D  *"
     print '(a)' ,"| (_| | |_| | (_| | (__| | | | (_) |      *"
     print '(a)' ," \__, |\__,_|\__,_|\___|_| |_|\___/       *"
     print '(a)' ," |___/                                    *"
  endif
#ifdef MPIP
  if(rank.eq.master) then
     print '(a,i3,a)','*    running with mpi in', np , ' processors    *'
     print '(a)' ,'*******************************************'
     print '(a)', 'Calculating Column Density'
  end if
  call mpi_cart_create(mpi_comm_world, ndim, dims, period, 1            &
       , comm3d, err)
  call mpi_comm_rank(comm3d, rank, err)
  call mpi_cart_coords(comm3d, rank, ndim, coords, err)
  print '(a,i3,a,3i4)', 'processor ', rank                              &
       ,' ready w/coords',coords(0),coords(1),coords(2)   
  call mpi_cart_shift(comm3d, 0, 1, left  , right, err)
  call mpi_cart_shift(comm3d, 1, 1, bottom, top  , err)
  call mpi_cart_shift(comm3d, 2, 1, out   , in   , err)
  call mpi_barrier(mpi_comm_world, err)   
  !
#else
  print '(a)' ,'*******************************************'
  print '(a)' ,'*     running on a single processor       *'
  print '(a)' ,'*******************************************'
  print '(a)', 'Calculating Column Density'
#endif

!   grid spacing
  dx=xmax/nxtot
  dy=ymax/nytot
  dz=zmax/nztot

!   allocate big arrays in memory
allocate( u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) )

end subroutine init_coldens

!=======================================================================

!> @brief reads data from file
!> @details reads data from file
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserved variables
!> @param integer [in] itprint : number of output
!> @param string [in] filepath : path where the output is

subroutine read_data(u,itprint,filepath)

  use parameters, only : &
                       np, neq, nxmin, nxmax, nymin, nymax, nzmin, nzmax
  use globals, only : rank, comm3d
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  integer, intent(in) :: itprint
  character (len=128), intent(in) :: filepath
  integer :: unitin, ip, err
  character (len=128) file1
  character           :: byte_read
  character, parameter  :: lf = char(10) 
  integer :: nxp, nyp, nzp, x0p, y0p, z0p, mpi_xp, mpi_yp, mpi_zp,neqp, neqdynp, nghostp
  real :: dxp, dyp, dzp, scal(3), cvp
  
  take_turns : do ip=0,np-1
    if (rank == ip) then


#ifdef MPIP
      write(file1,'(a,i3.3,a,i3.3,a)')  &
            trim(filepath)//'/BIN/points',rank,'.',itprint,'.bin'
      unitin=rank+10
#else
      write(file1,'(a,i3.3,a)')         &
            trim(filepath)//'/BIN/points',itprint,'.bin'
      unitin=10
#endif
      open(unit=unitin,file=file1,status='old', access='stream', &
            convert='LITTLE_ENDIAN')
   
      !   discard the ascii header
      do while (byte_read /= achar(255) )
        read(unitin) byte_read
        !print*, byte_read
      end do
      !  read bin header, sanity check to do
      read(unitin) byte_read
      read(unitin) byte_read
      read(unitin) nxp, nyp, nzp
      read(unitin) dxp, dyp, dzp
      read(unitin) x0p, y0p, z0p
      read(unitin) mpi_xp, mpi_yp, mpi_zp
      read(unitin) neqp, neqdynp
      read(unitin) nghostp
      read(unitin) scal(1:3)
      read(unitin) cvp
      read(unitin) u(:,:,:,:)
      close(unitin)

      print'(i3,a,a)',rank,' read: ',trim(file1)
        
#ifdef MPIP
      call mpi_barrier(comm3d,err)
#endif

    end if
  end do take_turns

end subroutine read_data

!=======================================================================

!> @brief gets position of a cell
!> @details Returns the position and spherical radius calculated with
!! respect to  the center of the grid
!> @param integer [in] i : cell index in the x direction
!> @param integer [in] j : cell index in the y direction
!> @param integer [in] k : cell index in the z direction
!> @param real    [in] x : x position in the grid 
!> @param real    [in] y : y position in the grid 
!> @param real    [in] z : z position in the grid 

  subroutine getXYZ(i,j,k,x,y,z)

    use globals,    only : dx, dy, dz, coords
    use parameters, only : nx, ny, nz, nxtot, nytot, nztot
    implicit none
    integer, intent(in)  :: i, j, k
    real,    intent(out) :: x, y, z
 
    x=(float(i+coords(0)*nx-nxtot/2)+0.5)*dx
    y=(float(j+coords(1)*ny-nytot/2)+0.5)*dy
    z=(float(k+coords(2)*nz-nztot/2)+0.5)*dz
        
  end subroutine getXYZ

!=======================================================================

!> @brief Rotation around the X axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid

subroutine rotation_x(theta,x,y,z,xn,yn,zn)

   ! rotation around the x axis by an angle theta

   implicit none
   real, intent(in ) :: theta, x, y, z
   reaL, intent(out) :: xn, yn, zn
   xn =   x
   yn =   y*cos(theta) - z*sin(theta)
   zn =   y*sin(theta) + z*cos(theta)
 end subroutine rotation_x

!=======================================================================

!> @brief Rotation around the Y axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid

 subroutine rotation_y(theta,x,y,z,xn,yn,zn)

   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) + z*sin(theta)
   yn =   y
   zn = - x*sin(theta) + z*cos(theta)
 end subroutine rotation_y

!=======================================================================

!> @brief Rotation around the Z axis
!> @details Does a rotation around the x axis
!> @param real [in], theta : Angle of rotation (in radians)
!> @param real [in], x : original x position in the grid
!> @param real [in], y : original y position in the grid
!> @param real [in], x : original z position in the grid
!> @param real [out], x : final x position in the grid
!> @param real [out], y : final y position in the grid
!> @param real [out], x : final z position in the grid

 subroutine rotation_z(theta,x,y,z,xn,yn,zn)

   implicit none
   real, intent(in ) :: theta, x, y, z
   real, intent(out) :: xn, yn, zn
   xn =   x*cos(theta) - y*sin(theta)
   yn =   x*sin(theta) + y*cos(theta)
   zn =   z
 end subroutine rotation_z

!=======================================================================

!> @brief Fill target map
!> @details Fills the target map of one MPI block
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param real [in] u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax) : 
!! conserved variables
!> @param real [out] map(nxmap,mymap) : Target map
!> @param real [in] dxT : target pixel width
!> @param real [in] dyT : target pixel height
!> @param real [in] thetax : Rotation around X
!> @param real [in] thetay : Rotation around Y
!> @param real [in] thetaz : Rotation around Z


subroutine fill_map(nxmap,nymap,u,map,dxT,dyT,theta_x,theta_y,theta_z)

   use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         neq, nx, ny, nz, vsc2, rsc,nztot, neqdyn, rhosc
  use globals, only : dz
  use hydro_core, only : u2prim

  implicit none

  integer, intent(in) :: nxmap,nymap
  real, intent(in) :: u(neq,nxmin:nxmax,nymin:nymax, nzmin:nzmax)
  real , intent(in) :: dxT, dyT, theta_x, theta_y, theta_z
  real, intent(out) :: map(nxmap,nymap)
  integer :: i,j,k, iobs, jobs
  real :: x,y,z,xn,yn,zn
  
  do k=1,nz
     do j=1,ny
        do i=1,nx

          !  obtain original position
          call getXYZ(i,j,k, x,y,z)
          
          !  do the rotation of the coordinates
          call rotation_x(theta_x,x,y,z,xn,yn,zn)
          call rotation_y(theta_y,xn,yn,zn,x,y,z)
          call rotation_z(theta_z,x,y,z,xn,yn,zn)
          ! This is the position on the target (centered)
          ! Integration is along Z
          iobs=xn/dxT + nxmap/2
          jobs=yn/dyT + nymap/2

          !  make sure the result lies in the map bounds
          if( (iobs >=1    ).and.(jobs >=1    ).and. &
              (iobs <=nxmap).and.(jobs <=nymap) ) then

            if (u(7,i,j,k) <= 0. ) then
              map(iobs, jobs) = map(iobs, jobs) + u(1,i,j,k)*rhosc*dz*rsc
            end if
          
          end if
      end do
    end do
  end do

end subroutine fill_map

!=======================================================================
!> @brief Writes header 
!> @details Writes header for binary input 
!> @param integer [in] unit : number of logical unit

subroutine write_header(unit, nx, ny)
  use globals, only : dx, dy, dz
  use parameters,  only : rsc, vsc, rhosc, cv
  implicit none
  integer, intent(in) :: unit, nx, ny
  character, parameter  :: lf = char(10) 
  character (len=128) :: cbuffer

  !  Write ASCII header
  write(unit) "**************** Output for Guacho v1.1****************",lf

  write(cbuffer,'("Dimensions    : ", i0,x,i0,x,i0)') NX, NY, 1
  write(unit) trim(cbuffer), lf
  
  write(cbuffer,'("Spacings      : ", 3(es10.3))') dX, dY, dZ
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Block Origin, cells    : ", i0,x,i0,x,i0)') 0, 0, 0
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("MPI blocks (X, Y, Z)   : ", i0,x,i0,x,i0)') 1, 1, 1
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Number of Equations/dynamical ones  ", i0,"/",i0)') 1, 1
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Number of Ghost Cells  ", i0)') 0
  write(unit) trim(cbuffer), lf

  write(cbuffer,'("Scalings ", a )')
  write(unit) trim(cbuffer), lf
 
  write(cbuffer, '("r_sc: ",es10.3," v_sc: ",es10.3," rho_sc: ",es10.3)') &
    rsc, vsc, rhosc
  write(unit) trim(cbuffer), lf

   write(cbuffer, '("Specfic heat at constant volume Cv: ",f7.2)') cv
  write(unit) trim(cbuffer), lf

#ifdef DOUBLEP
  write(unit) "Double precision 8 byte floats",lf
#else
  write(unit) "Double precision 4 byte floats",lf
#endif
 
  write(unit) "*******************************************************",lf
  write(unit) achar(255),lf
#ifdef DOUBLEP
  write(unit) 'd'
#else
  write(unit) 'f'
#endif
  write(unit) nx, ny, 1
  write(unit) dx, dy, dz
  write(unit) 0, 0, 0
  write(unit) 1, 1, 1
  write(unit) 1, 1
  write(unit) 0
  write(unit) rsc, vsc, rhosc
  write(unit) cv

end subroutine write_header


!=======================================================================

!> @brief Writes projection to file
!> @details Writes projection to file
!> @param integer [in] itprint : number of output
!> @param string [in] fileout : file where to write
!> @param integer [in] nxmap : Number of X cells in target
!> @param integer [in] nymap : Number of Y cells in target
!> @param real [in] map(nxmap,mymap) : Target map

subroutine  write_map(fileout,nxmap,nymap,map)

  implicit none
  integer, intent(in) :: nxmap, nymap
  character (len=128), intent(in) :: fileout
  real, intent(in) :: map(nxmap,nymap)
  integer ::  unitout

  unitout = 11
  open(unit=unitout,file=trim(fileout),status='unknown',access='stream', &
       convert='LITTLE_ENDIAN')

  call write_header(unitout,nxmap, nymap)

  write (unitout) map(:,:)
  close(unitout)
  
  print'(a,a)'," wrote file:",trim(fileout)
  
end subroutine write_map

!=======================================================================

end module coldens_utilities

!=======================================================================

!> @brief Computes the H-alpha emission
!> @details Computes the H-alpha apbsorption
!! @n It rotates the data along each of the coordinates axis
!! by an amount @f$ \theta_x, \theta_y, \theta_z @f$, and  projectcs the
!! map along the the LOS, which is taken to be the Z axis

program coldens

  use constants, only : pi
  use parameters, only : xmax,master, mpi_real_kind
  use globals, only : u, rank, comm3d
  use coldens_utilities
  implicit none
#ifdef MPIP
  include "mpif.h"
#endif
  character (len=128) :: pathin, fileout
  integer :: err
  integer :: itprint
  real    :: theta_x, theta_y, theta_z
  integer, parameter :: nxmap=192, nymap=512
  real :: dxT, dyT
  real :: map(nxmap, nymap), map1(nxmap,nymap)

  !  Target pixel size, relative to the simulation
  dxT= xmax/float(nxmap)
  dyT= dxT

  ! initializes program (uses the parameters.f90 form the rest of the code)
  call init_coldens()

  if (rank == master) then

    !  reads data form user
    print*,'Enter the input path (the path where the BIN directory is):  '
    read(5,*) pathin

    print*, 'Enter the number of output'
    read(5,*) itprint

    print*,' Enter the rotation angles theta_x, theta_y, theta_z (in degrees)'
    read(5,*) theta_x, theta_y, theta_z
    theta_x = theta_x * pi/180.
    theta_y = theta_y * pi/180.
    theta_z = theta_z * pi/180.

    print*,'Enter the output file: '
    read(5,*) fileout

  endif

  !  broadcasts user input to the rest of the processors
#ifdef MPIP
  call mpi_bcast(pathin ,128,mpi_character,0,mpi_comm_world,err)
  
  call mpi_bcast(itprint,  1,mpi_integer  ,0,mpi_comm_world,err)

  call mpi_bcast(theta_x,  1,mpi_real_kind,0,mpi_comm_world,err)
  call mpi_bcast(theta_y,  1,mpi_real_kind,0,mpi_comm_world,err)
  call mpi_bcast(theta_z,  1,mpi_real_kind,0,mpi_comm_world,err)

  call mpi_bcast(fileout,128,mpi_character,0,mpi_comm_world,err)

#endif

  !  read u from file
  call read_data(u,itprint,pathin)
  
  !  resets map
  map(:,:)=0.
  map1(:,:)=0.
  !
  if (rank == master) then
     print'(a)', 'Calculating projection with angles of rotaton'
     print'(f6.2,a,f6.2,a,f6.2,a)',theta_x*180./pi,'° around X, ' &
                                  ,theta_y*180./pi,'° around Y, '&
                                  ,theta_z*180./pi,'° around Z, '
  end if
  
  !  add info to the map
  call fill_map(nxmap,nymap,u,map,dxT,dyT, theta_x, theta_y, theta_z)
  !  sum all the partial sums
  call mpi_reduce(map,map1,nxmap*nymap, mpi_real_kind, mpi_sum, master, comm3d, err)
  
  !  write result
  if (rank == master) then 
    call write_map(fileout,nxmap,nymap,map1)
  end if

  if (rank == master) print*, 'my work here is done, have a  nice day'
#ifdef MPIP
  call mpi_finalize(err)
#endif
  !
  stop

end program coldens

!=======================================================================


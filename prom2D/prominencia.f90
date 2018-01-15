!=======================================================================
! Copyright (c) 2014 A. Esquivel, M. Schneiter, C. Villareal D'Angelo
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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.
!=======================================================================
module prominencia

  use parameters
  implicit none
  real :: RSW     !< Stellar radius
  real :: TSW     !< Stellar wind temperature
  real :: VSW     !< Stellar wind velocity
  real :: dsw     !< Stellar Wind Density
#if defined(PMHD) || defined(MHD)
  real :: b0                     !< Magnetic Field
#endif
  real :: MassS   !< Mass of the Star

contains

!=======================================================================

!> @brief Module initialization
!> @details Here the parameters of the Star are initialized, and scaled 
!! to code units

subroutine init_wind()
  
  use constants, only : Msun,Rsun, Yr, Ggrav, pi, au, day
  use globals, only :rank, dx
  implicit none
  real :: amdot  ! mdot_star (MSUN/yr)
#if defined(PMHD) || defined(MHD)
  real :: bsw, bsc
#endif

  !----------------STAR PARAMETERS ------------------
  MassS = msun
  AMDOT = 2.E-14*msun/yr              ! Stellar Mass Loss rate (g s^-1)

  TSW   = 3.0E6     !************      ! Stellar temperature (K)
  Vesc  = sqrt(2.*Ggrav*MassS/RSW)  !*******      !escape velocity (cm/s)
  a     = sqrt(2*Kb*tsw/(0.5*amh))
  

  RSW   = Rsun  !********************

  vsw   = a*(0.5*vesc/a)**2*exp(-0.5*(vesc/a)**2+3./2.)   !*************       ! Stellar wind velocity (cm/s)
  dsw   =((AMDOT/RSW)/(4*pi*RSW*VSW))   ! Stellar density @RS (g cm^-3)

  bsw    =1.0e3                         ! Stellar magnetic field (g)

  ! change to code units
  dsw=dsw/rhosc
  vsw=vsw/sqrt(vsc2)
  Tsw=Tsw/Tempsc
  Rsw=Rsw/rsc

#if defined(PMHD) || defined(MHD)
  bsw=bsw/bsc 
#endif

end subroutine init_wind

!=======================================================================

!> @brief Inject sources of wind
!> @details Imposes the sources of wond from the star and planet
!> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
!! conserver variables
!> @param real [time] time : current integration timr
  !--------------------------------------------------------------------

subroutine impose_wind(u,time)
  
  use constants, only : pi
  use globals, only : coords, dx, dy, dz, rank
  implicit none
  real, intent(out) :: u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax)
  real, intent (in) :: time
  real :: x, y, z, xpl, ypl, zpl
  real :: velx, vely, velz, rads, dens, radp, phi
  real :: vxorb, vyorb, vzorb, omega ! cpi
  integer ::  i,j,k


  do i=nxmin,nxmax
     do j=nymin,nymax
        do k=nzmin,nzmax

           ! Position measured in the grid
           x=(float(i+coords(0)*nx - nxtot/2.)+0.5)*dx
           y=(float(j+coords(1)*ny - nytot/2.)+0.5)*dy
           z=(float(k+coords(2)*nz - nztot/2.)+0.5)*dz

           ! Distance from the centre of the star
           rads=sqrt(x**2+y**2+z**2)

           ! in the bottom of the mesh we insert the wind 
           if( (y <= 0.2*Rsun) .and. (y>= Rsun) ) then

              VelX=0.
              VelY=VSW
              VelZ=0.
              DENS=DSW

              !   total density and momenta
              u(1,i,j,k) = dens
              u(2,i,j,k) = dens*velx
              u(3,i,j,k) = dens*vely
              u(4,i,j,k) = dens*velz

              !   magnetic field
#ifdef BFIELD 

              cpi=b0*(rsw/(rads+1.e-30))**3/(2.*(rads+1.e-30)**2)
              u(6,i,j,k) =  3.*z*x*cpi
              u(7,i,j,k) =  3.*y*z*cpi
              u(8,i,j,k) = (3.*z**2-rads**2)*cpi
#endif

	      if (mhd) then
#ifdef BFIELD
              ! total energy
              u(5,i,j,k)=0.5*dens*vsw**2         &
                   + cv*dens*2.*1.9999*Tsw       & 
                   + 0.5*(u(6,i,j,k)**2+u(7,i,j,k)**2+u(8,i,j,k)**2)
#endif
              endif
              
              if (passives) then
#ifdef PASSIVES
              !  density of neutrals
              u(neqdyn+1,i,j,k)= 0.0001*dens
              !   passive scalar (h-hot, c-cold, i-ionized, n-neutral)
              u(neqdyn+2,i,j,k)= dens   ! passive scalar
	      endif              
#endif                
           end if
              
        end do
     end do
  end do

end subroutine impose_wind

!=======================================================================
subroutine get_parker(r,Vrparker)

use constants, only: Ggrav,kB,amh

implicit none

real, intent(in)::r
real, intent(out)::vrparker
real:: vesc
real:: csS
real:: lambda
real:: Vr,Vr0, dVr,Rc,LHS, RHS, LHSold

!unscale the variables to work on physicals units
    vesc   = (2.*Ggrav*MassS/(Rsw*rsc))**0.5
    csS    = (kB*Tsw*Tempsc/(0.5*amh))**0.5
    lambda = 0.5*(vesc/csS)**2
! get sonic point radius
    Rc = Ggrav * MassS / (2.0 * csS**2.0)

    Vr0 = csS
    dVr = csS/10
    if (r < Rc) then
        dVr = - dVr
    end if

    LHS = Vr0 * exp(-Vr0**2.0 / (2.0 * csS**2.0) )
    RHS = csS * (Rc/r)**2.0 * exp(-2.0 * Rc/r + 3.0/2.0)

    Vr = Vr0
    do while (abs(LHS/RHS-1.0) > 1.e-8)
            ! save old LHS 
            LHSold = LHS
            !update Vr
            Vr = Vr + dVr
            ! calculate new LHS 
            LHS = Vr * exp(-Vr**2.0/(2.0 *csS**2.0) )
            ! see if crossed solution
            if((LHS/RHS-1.0)*(LHSold/RHS-1.0)<0)then
                dVr=-dVr/2.0
            end if
    end do
    Vrparker=Vr
end subroutine


!=======================================================================





end module prominencia

!=======================================================================

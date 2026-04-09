!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014

!! Description:
!!   calculates dissociation-/ionisation-/heating- rates for input into chemistry solver
!!
!! Input: 
!!   dtin   : timestep 
!!	 blockID: index of the block
!!   ind    : index of the zone inside the block
!!   dr     : ray segment length
!!   Nion   : number of ionizing photons on the ray
!!   Eion   : energy per photon
!!   Vpix   : volume of the zone

!! TODO implement proper stopping criterion based on ionization rate and accumulated column Av, density
!! TODO continuous shielding function from Draine & Bertoldi 1996, 
!!      as otherwise not fully photon conservative, jump in shielding leads to different behaviour
!!      depending if the column is small enough or not to see the jump, only issue for very high res.


!! IMPORTANT for no H2 ionisation set Nih2 AND sH2i to 0

!#define DEBUG_RAY
!#define SANITY
#define PROC 0

subroutine pt_solveZone(dtin, blockID, ind, dr, Eion, sH, Nion, Vpix, stopp, dirx, diry, dirz)
  
  ! access to state variables
  use Grid_interface, ONLY : Grid_getPointData, Grid_putPointData
  ! constants
  use Particles_rayData, ONLY :  ph_radPressure, speedoflight

  use rt_data, ONLY : rt_abar, rt_protonMass

  !  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  
#include "Flash.h"
#include "constants.h"

  
  ! zone position
  integer, intent(in) :: blockID
  logical, intent(out) :: stopp
  integer, dimension(NDIM),intent(in) :: ind
  real, intent(in) :: dr, dtin, Vpix

  real, intent(inout) :: Eion, sH
  real, intent(inout) :: Nion
  real, intent(in) :: dirx, diry, dirz

  real :: tmp, xH0, dens, numdens, xHp
  real :: phih, hvphih

  real :: DtauHI, DNionHI
! radiation pressure
  real :: mh, velx, vely, velz, HionMom
  real :: FullEion
  logical :: stoppH

! those are always values for the zone centers
! get zone data, rho, neutral hydrogen, temperature
  call Grid_getPointData(blockID, CENTER, DENS_VAR, INTERIOR, ind, dens)

  call Grid_getPointData(blockID, CENTER, IHA_SPEC, INTERIOR, ind, xH0)
  call Grid_getPointData(blockID, CENTER, IHP_SPEC, INTERIOR, ind, xHp)

  if(xH0 .eq. 0.0 .or. xHp .eq. 1.0) then
! need something to calculate otherwise enjoy some NaNs
   xH0 = 1e-10
  endif

  HionMom  = 0.0
  stopp   = .false.

! number density per H atom
  numdens = dens/(rt_abar*rt_protonMass)

! photo-heating rate 13.6 to 15.2
  call Grid_getPointData(blockID, CENTER, PHHE_VAR, INTERIOR, ind, hvphih)
! ionization rate 13.6 to 15.2 eV
  call Grid_getPointData(blockID, CENTER, PHIO_VAR, INTERIOR, ind, phih)

! weight cross section by median 
  DtauHI = numdens * dr * sH * xH0 

! number of photons that ionize
  DNionHI = Nion*(1d0-dexp(-DtauHI))

! Here we must add the ionization energy back to the energy
! of the photon before calculating momentum.
! Currently this is just hacking it in, but it should
! be a runtime parameter passed through the par file. -JW
! Note here I'm assuming we are ionizing hydrogen. -JW

!13.6 eV*1 eV in ergs + Thermal energy imparted to electon. - JW
  FullEion = 13.6*1.6022d-12 + Eion 

! m_{ph} = N photons * energy of photons / c
  HionMom = DNionHI*FullEion/speedoflight

! ionisation rate
! number of ionizing photons per timestep divided by number of neutral hydrogen atoms 
! N_ion(dt)/N_neutral
  tmp  = DNionHI/(numdens*xH0*Vpix*dtin)
  phih = phih + tmp

! photoheating rate, in units of erg/(cm^3 s) 
  hvphih = hvphih + tmp*Eion*numdens*xH0
  
  call Grid_putPointData(blockID, CENTER, PHIO_VAR, INTERIOR, ind, phih)
  call Grid_putPointData(blockID, CENTER, PHHE_VAR, INTERIOR, ind, hvphih)
 
  Nion = Nion - DNionHI

  if(Nion .le. 1.0) then
    Nion   = 0.0 
    stopp = .true.
  endif

! Radiation pressure 
  if(ph_radPressure) then
    mh = HionMom/(dens*Vpix)

! put back to data structure
    call Grid_getPointData(blockID, CENTER, VELX_VAR, INTERIOR, ind, velx)
    call Grid_getPointData(blockID, CENTER, VELY_VAR, INTERIOR, ind, vely)
    call Grid_getPointData(blockID, CENTER, VELZ_VAR, INTERIOR, ind, velz)

! add momentum
    velx = velx + mh*dirx
    vely = vely + mh*diry
    velz = velz + mh*dirz

    call Grid_putPointData(blockID, CENTER, VELX_VAR, INTERIOR, ind, velx)
    call Grid_putPointData(blockID, CENTER, VELY_VAR, INTERIOR, ind, vely)
    call Grid_putPointData(blockID, CENTER, VELZ_VAR, INTERIOR, ind, velz)
  endif

  return
end subroutine pt_solveZone

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
!#define ONE_CELL_TESTING
#define PROC 0

subroutine pt_solveZone(dtin, blockID, ind, dr, Eion, sH, Nion, Vpix, &
                        stopp, dirx, diry, dirz, ph_type, zone_size)
  
  ! access to state variables
  use Grid_interface, ONLY : Grid_getPointData, Grid_putPointData
  ! constants
  use Particles_rayData, ONLY :  ph_radPressure, speedoflight, &
              ion_photon, pe_photon, early_term_FUV, ph_EUVonDust

  use rt_data, ONLY : rt_abar, rt_protonMass, sigDust, dust_gas_ratio, rt_dt ! global ray trace dt
  use Heat_data, ONLY : he_boltz, he_dust_sputter_temp

  !  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  
#include "Flash.h"
#include "constants.h"

  
  ! zone position
  integer, intent(in) :: blockID, ph_type
  logical, intent(out) :: stopp
  integer, dimension(NDIM),intent(in) :: ind
  real, intent(in) :: dr, dtin, Vpix, zone_size

  real, intent(inout) :: Eion, sH
  real, intent(inout) :: Nion
  real, intent(in) :: dirx, diry, dirz

  real :: tmp, xH0, dens, numdens, xHp
  real :: phih, hvphih, Flux, GFlux

  real :: DtauHI, DNionHI
  real :: DtauDust, DNdust ! Maximum num photons abs by H, rest go to dust.
! radiation pressure
  real :: velx, vely, velz, HionMom, DustMom, HionVel, DustVel, delVel
  real :: FullEion, NumH0
  logical :: stoppH, fully_ionized
  real, parameter   :: Newton = 6.6725985d-8, pi = 3.1415926535897932384d0, N_H0=1.87e21
  real    :: Av, f_ext, lam_j, temp, cellFlux, bgFlux
  real    :: DNambient, ambientFlux
  real    :: surface_correction, atan2xy, costheta2

! if dr = 0, we didn't move the ray so just return.

if (dr .le. 0.0d0) return

! those are always values for the zone centers
! get zone data, rho, neutral hydrogen, temperature
  call Grid_getPointData(blockID, CENTER, DENS_VAR, INTERIOR, ind, dens)
! Neutral fraction.
  call Grid_getPointData(blockID, CENTER, IHA_SPEC, INTERIOR, ind, xH0)
! Ionized fraction.
  call Grid_getPointData(blockID, CENTER, IHP_SPEC, INTERIOR, ind, xHp)
! ionization rate 13.6 to 15.2 eV
  call Grid_getPointData(blockID, CENTER, PHIO_VAR, INTERIOR, ind, phih)

  
  DNionHI = 0.0d0
  DNdust  = 0.0d0
  Flux    = 0.0d0
  stopp   = .false.
  HionMom = 0.0d0
  DustMom = 0.0d0
  HionVel = 0.0d0
  DustVel = 0.0d0
  delVel  = 0.0d0

! number density per H atom
  numdens = dens/(rt_abar*rt_protonMass)
! Number of total neutral hydrogen in this cell.
  NumH0 = xH0*Vpix*numdens
  
! You cannot photoionize more hydrogen than neutral hydrogen exists in the cell!
  !fully_ionized = (phih*rt_dt .ge. NumH0) .or. (xHp .ge. 1.0d0) .or. (xH0 .le. 0.0d0)
  fully_ionized = ((xHp .ge. 1.0d0) .or. (xH0 .le. 0.0d0))
! print*, "[pt_SolveZone]: photon_type =", ph_type, ion_photon, pe_photon

#ifdef ONE_CELL_TESTING
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: Begin phi*rt_dt =", phih*rt_dt
            write(*,*) "[pt_solveZone]: fully_ionized=", fully_ionized
#endif



!surface_correction = 1.
! First order correction to effective surfaces of cells -MW
costheta2 = dirz*dirz/(dirx*dirx + diry*diry + dirz*dirz)
atan2xy = atan2(dirx, diry)

surface_correction = sqrt(1.-costheta2) + sqrt(costheta2)
surface_correction = surface_correction*(abs(sin(atan2xy)) + abs(cos(atan2xy)))


!phi = atan2(diry, dirx)
!theta = acos(dirz/sqrt(dirx*dirx + diry*diry + dirz*dirz))

!surface_correction = abs(sin(phi)) + abs(cos(phi))
!surface_correction = surface_correction*(abs(sin(theta)) + abs(cos(theta)))

!if (phi < -pi/2.) then
!  surface_correction = -sin(phi) - cos(phi)
!else if (phi < 0.) then
!  surface_correction = -sin(phi) + cos(phi)
!else if (phi < pi/2.) then
!  surface_correction =  sin(phi) + cos(phi)
!else
!  surface_correction =  sin(phi) - cos(phi)
!end if
!if (theta < -pi/2.) then
!  surface_correction = surface_correction * (-sin(theta) - cos(theta))
!else if (theta < 0.) then
!  surface_correction = surface_correction * (-sin(theta) + cos(theta))
!else if (theta < pi/2.) then
!  surface_correction = surface_correction * ( sin(theta) + cos(theta))
!else
!  surface_correction = surface_correction * ( sin(theta) - cos(theta))
!end if


! Select the correct bin before calculating ionization and heating. - JW
if ( (ph_type == ion_photon) .and. (.not. fully_ionized) ) then

!!!!!!!!
!! NOTE NOTE NOTE Need to add in dust absorption (with dust cross section).
!! Then can add this heating to Gflux after.
!!!!!!!!

! photo-heating rate 13.6 to 15.2
  call Grid_getPointData(blockID, CENTER, PHHE_VAR, INTERIOR, ind, hvphih)

  if(xH0 .eq. 0.0d0 .or. xHp .eq. 1.0d0) then
! need something to calculate otherwise enjoy some NaNs
  print*, "[pt_solveZone]: We shouldn't ever be in here again."
   xH0 = 1e-50
  endif
  
!  write(*,'(A,9ES10.3)') "[pt_solveZone]: ephen, kion, xH0, xHp, Nion, &
!                            DNionHI, ndens, dr, sH =", hvphih, phih, xH0, xHp, Nion, &
!                                                   DNionHI, numdens, dr, sH
  
! weight cross section by median 
  DtauHI = numdens * dr * sH * xH0 

! number of photons that ionize
  DNionHI = Nion*(1d0-dexp(-DtauHI))

!if (DNionHI .gt. NumH0) then
!  write(*,*) "[pt_solveZone]: Hey, that's too many photons!"
!  write(*,'(A,ES13.3E3)') "[pt_solveZone]: DNionHI =", DNionHI
!  write(*,'(A,ES13.3E3)') "[pt_solveZone]: NumH0 =", NumH0
  
!  DNionHI = min(DNionHI, NumH0)
!end if

#ifdef ONE_CELL_TESTING
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: Begin DNionHI =", DNionHI
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: Begin ion =", xHp
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: Begin neutral =", xH0
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: Begin photo ion rate =", phih
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: Begin photo heating rate =", hvphih
            call flush(6)
        !stop
#endif
! ionisation rate
! number of ionizing photons per timestep divided by number of neutral hydrogen atoms 
! N_ion(dt)/N_neutral
  tmp  = DNionHI/(numdens*xH0*Vpix*dtin)
  phih = phih + tmp

! photoheating rate, in units of erg/(cm^3 s) 
  hvphih = hvphih + tmp*Eion*numdens*xH0
  
  if (isNaN(phih)) then
    write(*,'(A,8ES10.3)') "[pt_solveZone]: ephen, kion, xH0, xHp, Nion, &
                            DNionHI, ndens, dr =", hvphih, phih, xH0, xHp, Nion, &
                                                   DNionHI, numdens, dr
    stop
  end if
  
  call Grid_putPointData(blockID, CENTER, PHIO_VAR, INTERIOR, ind, phih)
  call Grid_putPointData(blockID, CENTER, PHHE_VAR, INTERIOR, ind, hvphih)
  
! Here we must add the ionization energy back to the energy
! of the photon before calculating momentum.
! Currently this is just hacking it in, but it should
! be a runtime parameter passed through the par file. -JW
! Note here I'm assuming we are ionizing hydrogen. -JW

! Also note that this full Eion is what should be passed on to
! interact with dust down below.

!13.6 eV*1 eV in ergs + Thermal energy imparted to electon. - JW
  FullEion = 13.6*1.6022d-12 + Eion 

! Get the flux in scratch variable.
!#ifdef UVFL_VAR
!  call Grid_getPointData(blockID, CENTER, UVFL_VAR, INTERIOR, ind, Flux)
!#endif

  Flux = phih*numdens*xH0*dr*FullEion
  !Flux = phih*numdens*xH0*FullEion * Vpix**(1./.3) / surface_correction
  
! Store the flux in a scratch variable to look at later in plot files.
#ifdef UVFL_VAR
  call Grid_getPointData(blockID, CENTER, UVFL_VAR, INTERIOR, ind, cellFlux)
  Flux = Flux + cellFlux
  call Grid_putPointData(blockID, CENTER, UVFL_VAR, INTERIOR, ind, Flux)

  ! Store unabsorbed (ambient) flux similarly. - MW
  Flux = (Nion-DNionHI)*dr/(dtin*Vpix) * FullEion
  !Flux = (Nion-DNionHI)/(dtin*Vpix**(2./3.)*surface_correction) * FullEion
  call Grid_getPointData(blockID, CENTER, AUVF_VAR, INTERIOR, ind, cellFlux)
  Flux = Flux + cellFlux
  call Grid_putPointData(blockID, CENTER, AUVF_VAR, INTERIOR, ind, Flux)
#endif
  Flux = 0.0d0

! m_{ph} = N photons * energy of photons / c
  HionMom = DNionHI*FullEion/speedoflight
  HionVel = HionMom/(dens*Vpix)
  
! Remove the photons that ionized hydrogen.
  Nion = Nion - DNionHI


else

  FullEion = Eion ! The energy is already all there.

end if  ! else if (ph_type == pe_photon) then

! Now we just let everything interact with dust. Since ionizing photons
! have a much larger cross section with hydrogen, this effect only becomes
! important for them in ionized gas. Therefore we just always do this
! after photons are absorbed by hydrogen.

! Note: If the dust is sputtered, then clearly nothing can be absorbed
! by it. So check if we are above the sputtering temperaute and if so
! the ray passes through without hitting any dust. - JW
call Grid_getPointData(blockID, CENTER, TEMP_VAR, INTERIOR, ind, temp)

if ((Nion .gt. 1.0d0) .and. (temp < he_dust_sputter_temp) .and. &
    (ph_type == pe_photon .or. ph_EUVonDust)) then

#ifdef ONE_CELL_TESTING
  print*, "Am I a UV photon?", (ph_type == ion_photon)
#endif

! I'm lazy, I didn't rename the Nion variable. Just note that for radiation
! bin for PE photons (5.6-13.6 eV) when you see Nion it means pe photons. - JW

  call Grid_getPointData(blockID, CENTER, PEFL_VAR, INTERIOR, ind, GFlux)

! weight cross section by median

! Note for photoelectric effect on dust we use the cross section of dust
! per H nucleon sigma=1e-21 [cm^2 / H] from Draine 2011.

! Note sigma = n_d / n_H * regular cross section dust [cm^2]
! and so it has the number fraction of dust grains already included, 
! and therefore doesn't need the metallicity in the equation. - JW
  DtauDust = numdens * dr * sigDust ! * 0.02

! number of photons that strike dust grains. - JW
  DNdust = Nion*(1d0-dexp(-DtauDust))

! the rest of the photons, which form the ambient FUV field. - MW
  DNambient = Nion - DNdust
  
! Convert this number into a flux by multiplying by the
! average energy of a photon in this bin. - JW

  Flux = DNdust*FullEion/dtin
  ambientFlux = DNambient*FullEion/dtin

! Now convert this to a Habing 1968 (G_conv = 1.6 x 10^-3 ergs cm^-2 s^-1) normalized flux as
! described in Baczynski 2015. Note we can probably do
! better if we used the actual cell entry angle instead of face area, but
! good enough for now.

! Also note, implict assumption here that cells are cubes! - JW
  Flux = Flux / (1.6d-3 * surface_correction*zone_size**2.0)
  !Flux = Flux / (1.6d-3 * Vpix/dr)
  GFlux = GFlux + Flux

  ambientFlux = ambientFlux / (1.6d-3 * surface_correction*zone_size**2.0)
  !ambientFlux = ambientFlux / (1.6d-3 * Vpix/dr)

  if (isNaN(GFlux)) then
    write(*,'(A,9ES10.3)') "[pt_solveZone]: ephen, kion, xH0, xHp, Nion, DNionHI, &
                            ndens, dr, Flux =", hvphih, phih, xH0, xHp, Nion, &
                            DNionHI, numdens, dr, Flux
    stop
  end if

! This flux is then stored in the cell to be used to calculate the proper
! heating from PE effect in combination with other heating and cooling
! sources in RadHeat.F90 - JW

  call Grid_putPointData(blockID, CENTER, PEFL_VAR, INTERIOR, ind, GFlux)

#ifdef FUFL_VAR
! If this is FUV flux and the background FUV has become the dominant flux,
! terminate the ray (otherwise the code will follow these rays literally for
! parsecs while they are orders of magnitude smaller than the background FUV). - JW
  if (early_term_FUV) then
    if (ph_type == pe_photon) then
    
        lam_j  = 1.0/(dens/numdens) * sqrt(5./3.*pi*he_boltz*temp / numdens / Newton)
        Av     = lam_j * numdens * dust_gas_ratio / N_H0 ! Using visual extinction estimation.
        f_ext  = exp(-3.5*Av) ! Fraction of background FUV that makes it through to here.
        !Av     = lam_j * numdens * sigDust      ! Using actual dust cross section.
        !f_ext  = exp(-Av) ! Fraction of background FUV that makes it through to here.
        bgFlux = 1.69*f_ext! Draine field is 1.69*G_0
    
        ! If flux on this cell and flux still in the ray are less than the background
        ! Draine field, terminate this ray.
        if (Flux < bgFlux .and. ambientFlux < bgFlux) then 
            !print*, "bgFlux =", bgFlux, "Flux =", Flux
            Nion = 0.0
            stopp = .true.
        end if
    end if
  end if
  
! Store the flux in a scratch variable to look at later in plot files.
  if (ph_type == pe_photon) then
    call Grid_getPointData(blockID, CENTER, FUFL_VAR, INTERIOR, ind, cellFlux)
    Flux = Flux*1.6d-3 + cellFlux ! Convert back to ergs cm^-2 s^-1
    call Grid_putPointData(blockID, CENTER, FUFL_VAR, INTERIOR, ind, Flux)

    ! Unabsorbed, ambient flux, stored similarly.
    call Grid_getPointData(blockID, CENTER, AFUF_VAR, INTERIOR, ind, cellFlux)
    ambientFlux = ambientFlux*1.6d-3 + cellFlux
    call Grid_putPointData(blockID, CENTER, AFUF_VAR, INTERIOR, ind, ambientFlux)
  end if
#endif

! Add momentum imparted on the dust by photons. - JW
! m_{ph} = N photons * energy of photons / c
  DustMom = DNdust*FullEion/speedoflight
  DustVel = DustMom/(dens*Vpix)

  Nion = Nion - DNdust

  
end if


  if(Nion .le. 1.0) then
    Nion   = 0.0 
    stopp = .true.
  endif

! Radiation pressure 
  if(ph_radPressure) then

    delVel = HionVel + DustVel
    
!    write(*,'(A,ES10.3)') "HionMom=", HionMom
!    write(*,'(A,ES10.3)') "DustMom", DustMom
!    write(*,'(A,ES10.3)') "HionVel=", HionVel
!    write(*,'(A,ES10.3)') "DustVel=", DustVel
!    write(*,'(A,ES10.3)') "delVel=", delVel
!    call flush(6)
    !stop

! put back to data structure
    call Grid_getPointData(blockID, CENTER, VELX_VAR, INTERIOR, ind, velx)
    call Grid_getPointData(blockID, CENTER, VELY_VAR, INTERIOR, ind, vely)
    call Grid_getPointData(blockID, CENTER, VELZ_VAR, INTERIOR, ind, velz)

! add momentum
    velx = velx + delVel*dirx
    vely = vely + delVel*diry
    velz = velz + delVel*dirz

    call Grid_putPointData(blockID, CENTER, VELX_VAR, INTERIOR, ind, velx)
    call Grid_putPointData(blockID, CENTER, VELY_VAR, INTERIOR, ind, vely)
    call Grid_putPointData(blockID, CENTER, VELZ_VAR, INTERIOR, ind, velz)
  endif

#ifdef ONE_CELL_TESTING
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: End Flux =", Flux
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: End temp =", temp
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: End Nion =", Nion
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: End photo ion rate =", phih
            write(*,'(A,ES13.3E3)') "[pt_solveZone]: End photo heating rate =", hvphih
            call flush(6)
        !stop
#endif

  return
end subroutine pt_solveZone

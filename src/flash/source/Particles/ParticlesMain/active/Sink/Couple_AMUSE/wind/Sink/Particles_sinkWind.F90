!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! subroutine Particles_sinkWind()
!!
!! Calls the wind injection subroutine over all the sink
!! particles. 
!!
!! Still to do: Add runtime parameter for the min mass here!
!!
!! Joshua Wall 5/26/2016
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define DEBUG

subroutine Particles_sinkWind(dt)

use Particles_interface, ONLY: Particles_wind
use Grid_interface, only : Grid_fillGuardCells
use Particles_sinkData, ONLY: localnpf, particles_global, min_wind_dt, &
                              localnp, particles_local
use Driver_data, ONLY: dr_globalMe, dr_simTime

use pt_sinkinterface, ONLY: pt_sinkGatherGlobal
use Particles_sort
use RuntimeParameters_interface, ONLY: RuntimeParameters_get

#include "Flash.h"
#include "constants.h"

real, intent(in)        :: dt

real                    :: dmdt, mass, max_dt, save_dt, inj_mass
real                    :: xloc, yloc, zloc, inj_vel_mag, twind, bgdy
real, save              :: min_wind_mass
real, allocatable       :: mass_sorted(:)
integer, allocatable    :: QSindex(:)
integer                 :: i
logical, save           :: first_call = .true.

if (first_call) then
  call RuntimeParameters_get("min_wind_mass", min_wind_mass)
  first_call = .false.
end if

save_dt = 1e99

#ifdef DEBUG
  print*, "In Particles_sinkWind."
#endif

call pt_sinkGatherGlobal()

!!! MUY IMPORTANTE !!!
!!! We must sort the particles_global array the same
!!! on all procs for this to work correctly. - JW

allocate(mass_sorted(localnpf))
allocate(QSindex(localnpf))

do i=1, localnpf
  mass_sorted(i) = particles_global(MASS_PART_PROP,i)
end do

call sortInd(mass_sorted, QSindex)

#ifdef DEBUG
do i=1, localnpf
  print*, i, " on proc ", dr_globalMe, "mass =", &
          particles_global(MASS_PART_PROP,QSindex(i)) / 1.989e33, "Msun."
end do
#endif

do i=1, localnpf
  
  mass = particles_global(MASS_PART_PROP, QSindex(i))
  dmdt = particles_global(DMDT_PART_PROP, QSindex(i))
  ! Get background gas density. - JW
  bgdy = particles_global(BGDY_PART_PROP, QSindex(i))

    if ((mass .ge. min_wind_mass) .and. dmdt .gt. 0.0) then

      max_dt  = dt
    ! Fit of wind velocity at r = inf vs. star mass, from Dale et al. 2013
      !inj_vel_mag = 1.019430 * ( mass - 3.579183e34 )**0.24 + 6e7
      !inj_vel_mag = 1.0d8 !2.4d8
      inj_vel_mag = particles_global(VELW_PART_PROP, QSindex(i)) 
      dmdt = particles_global(DMDT_PART_PROP, QSindex(i))
      xloc = particles_global(POSX_PART_PROP, QSindex(i))
      yloc = particles_global(POSY_PART_PROP, QSindex(i))
      zloc = particles_global(POSZ_PART_PROP, QSindex(i))
      twind = dr_simTime + dt - particles_global(CREATION_TIME_PART_PROP, QSindex(i)) 
      !call Particles_wind(mass, dmdt, xloc, yloc, zloc, max_dt)

      inj_mass = dt*dmdt
      !inj_mass = dt * 1d-6 * 1.989d33 / (60d0*60d0*24d0*365d0)

      ! Wolf-Rayet in Fryer 2006
      !inj_mass = dt * 3d-3 * 1.989d33 / (60d0*60d0*24d0*365.25d0)
      ! ~30 MSun from Weaver 1977 
      !inj_mass = dt * 1d-6 * 1.989d33 / (60d0*60d0*24d0*365.25d0)
      
      ! subroutine inject_direct may now modify mass and velocity (while conserving
      ! momentum), so these write statments may not be valid
      !if (dr_globalMe .eq. MASTER_PE) then
      !write(*,'(A,ES13.3e3)') "Inject mass =", inj_mass
      !write(*,'(A,ES13.3e3)') "Inject vel =", inj_vel_mag
      !end if
      
#ifdef DEBUG
        print*, "Mass of sink =", particles_global(MASS_PART_PROP, QSindex(i)) / 1.989e33, "Msun."
        print*, "Blk prop =", particles_global(BLK_PART_PROP,QSindex(i))
        print*, "DMDT =", dmdt / 1.989d33 * (60d0*60d0*24d0*365.25d0), "Msun / yr."
        print*, "Wind vel =", inj_vel_mag / 1e5, "km / s."
        print*, "[Particles_sinkWind]: background density =", bgdy, "g / cm^3."
#endif
      
      
        call inject_direct([xloc, yloc, zloc], inj_mass, inj_vel_mag,  twind, dt, bgdy)
        ! Remove mass argument which will no longer work with inject_direct. - SA 20240408
        ! However, the star particle inject_direct (which has been updated to account
        ! for and implement jets) has several arguments not used here.
        ! Adding them will require defining the angular momentum and jet_wind switch for the
        ! sink particles within this routine.
#ifdef DEBUG
        print*, "[Particles_sinkWind]: background density =", bgdy, "g / cm^3."
#endif

        particles_global(BGDY_PART_PROP, QSindex(i)) = bgdy

        if (max_dt .lt. save_dt) save_dt = max_dt

    end if

end do

deallocate(mass_sorted)
deallocate(QSindex)

do i=1, localnp

    particles_local(BGDY_PART_PROP,i) = particles_global(BGDY_PART_PROP,i)
    
end do

call Grid_fillGuardCells(CENTER, ALLDIR, doEos=.false.)

!if (save_dt .lt. dt) then 
!  min_wind_dt = save_dt
!else
!  min_wind_dt = 1d99
!end if

end subroutine Particles_sinkWind

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
use Particles_sinkData, ONLY: localnpf, particles_global, min_wind_dt
use Driver_data, ONLY: dr_globalMe, dr_simTime

#include "Flash.h"
#include "constants.h"

real, intent(in) :: dt

real    :: dmdt, mass, max_dt, save_dt, inj_mass
real    :: xloc, yloc, zloc, inj_vel_mag, twind
integer :: i

save_dt = 1e99

#ifdef DEBUG
  print*, "In Particles_sinkWind."
#endif

!!! MUY IMPORTANTE !!!
!!! We must sort the particles_global array the same
!!! on all procs for this to work correctly. - JW

do i=1, localnpf
  
  mass = particles_global(MASS_PART_PROP,i)

#ifdef DEBUG
  print*, "Mass of sink =", mass
  print*, "Blk prop =", particles_global(BLK_PART_PROP,i)
#endif

  if (mass .ge. 8.0d0*1.989d33) then

  max_dt  = dt
! Fit of wind velocity at r = inf vs. star mass, from Dale et al. 2013
  !inj_vel_mag = 1.019430 * ( mass - 3.579183e34 )**0.24 + 6e7
  inj_vel_mag = 1.0d8
  dmdt = particles_global(DMDT_PART_PROP, i)
  xloc = particles_global(POSX_PART_PROP, i)
  yloc = particles_global(POSY_PART_PROP, i)
  zloc = particles_global(POSZ_PART_PROP, i)
  twind = dr_simTime + dt - particles_global(CREATION_PART_PROP, i) 
  !call Particles_wind(mass, dmdt, xloc, yloc, zloc, max_dt)

  !inj_mass = dt*dmdt
  inj_mass = dt * 1d-6 * 1.989d33 / (60d0*60d0*24d0*365d0)
  
  if (dr_globalMe .eq. MASTER_PE) then
  write(*,'(A,ES13.3e3)') "Inject mass =", inj_mass
  write(*,'(A,ES13.3e3)') "Inject vel =", inj_vel_mag
  end if

  call inject_direct([xloc, yloc, zloc], inj_mass, inj_vel_mag, twind, dt)

  if (max_dt .lt. save_dt) save_dt = max_dt

  end if

end do

if (save_dt .lt. dt) then 
  min_wind_dt = save_dt
else
  min_wind_dt = 1d99
end if

end subroutine Particles_sinkWind

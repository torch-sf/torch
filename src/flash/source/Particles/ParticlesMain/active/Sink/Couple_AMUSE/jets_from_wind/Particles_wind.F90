!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! subroutine Particles_wind
!!!
!!!
!!! Injects a stellar wind from a particle or sink particle given that
!!! particle's dM/dt and wind terminal velocity v_wind.
!!!
!!! J. Wall Drexel University 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define debug
#define debug2
#define debug_jets

subroutine Particles_wind(dt)

use Particles_data, only : particles, pt_typeInfo, pt_numLocal
use Particles_windData
use Particles_interface, only : Particles_getGlobalNum

!use Particles_windData

use Driver_interface, only : Driver_abortFlash
use Driver_data, only : dr_globalComm, dr_globalNumProcs, dr_globalMe, dr_simTime

use Grid_interface, only : Grid_fillGuardCells

use RuntimeParameters_interface, ONLY: RuntimeParameters_get

implicit none

real, intent(in)       :: dt

! Injection info for winds.
! Global vars
real*8, allocatable :: x(:), y(:), z(:), &
               j_x(:), j_y(:), j_z(:), &
               dmdt(:), v_wind(:), c_time(:), bgdy(:)
! Added ang momentum -SA 20230720

integer    :: w_num
! Local vars
real, allocatable      :: locx(:), locy(:), locz(:), locdmdt(:), locv_wind(:), &
			  locc_time(:), locbgdy(:)

real, allocatable      :: angmom_x(:), angmom_y(:), angmom_z(:) ! Added ang momentum -SA 20230720

real                   :: mass, twind, bgdy_old, time
			  
! For counting particles.
integer                :: p_begin, p_end, p_num, p_globalnum, w_numloc

! Storage of the local indices that are actual wind stars.
integer, allocatable   :: p_ind(:)

! Handling jet versus wind option -SA 20230726
logical :: jet_switch
integer, allocatable :: jet_wind(:), jw_switch(:)
integer, parameter :: jet_flag = 1 ! update -SA 20230802
integer, parameter :: wind_flag = 2

! For MPI comm
integer                :: num_array(dr_globalNumProcs), &
                          disp(dr_globalNumProcs), ierr, i, p, &
			  rank_minus_one

! Runtime parameters
character(len=8), save :: part_type
real, parameter :: yr = (60.0d0**2.0)*24.0d0*365.25d0
real, parameter :: solarMass = 1.989d33
logical, save          :: first_call = .true.

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"
#include "Particles.h"
#include "Eos.h"

if (first_call) then
  call RuntimeParameters_get("min_wind_mass", min_wind_mass)
  call RuntimeParameters_get("min_jet_mass", min_jet_mass)
  call RuntimeParameters_get("max_jet_mass", max_jet_mass)
  call RuntimeParameters_get("jet_time", jet_time)
  first_call = .false.
end if

! Access current time of the code (see also twind below):
time = dr_simTime + dt  !SA 20230728

min_wind_dt = 1d99


! Local number of massive/active particles.
  p_begin = pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE)
  p_num   = pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE)
  p_end   = p_num + p_begin - 1

! print*, "Number of particles (begin, num, end): ", p_begin, p_num, p_end, " -SA 202212"

allocate(p_ind(pt_numLocal))
p_ind = 0

call Particles_getGlobalNum(p_globalnum)
! Why are we getting the global number of particles here????
! The do loop that assigns values to these arrays uses the number on the processor....? -SA
! print*, "Global number of particles: ", p_globalnum
allocate(locx(p_globalnum), locy(p_globalnum), locz(p_globalnum), locc_time(p_globalnum))
allocate(angmom_x(p_globalnum), angmom_y(p_globalnum), angmom_z(p_globalnum))  !Added ang momentum -SA 20230720
allocate(locdmdt(p_globalnum), locv_wind(p_globalnum), locbgdy(p_globalnum))
allocate(jet_wind(p_globalnum)) ! Added 20230726 -SA

w_numloc  = 0
w_num     = 0
num_array = 0

locx = 0.0d0; locy=0.0d0; locz=0.0d0
angmom_x = 0.0d0; angmom_y = 0.0d0; angmom_z = 0.0d0 !Added ang momentum -SA 20230720
locdmdt = 0.0d0; locv_wind=0.0d0; locbgdy=0.0d0; locc_time= 0.0d0
jet_wind = 0 ! Added 20230726 -SA

#ifdef debug_jets
print*, "Before do loop - check arrays:", locx, angmom_x, jet_wind
print*, "Particles on this proces: ", p_num, "Proces: ", dr_globalMe
call flush()
#endif

do p = p_begin, p_end
#ifdef debug
  print*, "Particle mass =", particles(MASS_PART_PROP, p)
  print*, "Particle dmdt =", particles(DMDT_PART_PROP, p)
#endif 
  w_numloc = w_numloc + 1  !Moved to start of loop -SA 20230811

  ! Now testing for jet and wind conditions. SA 20230728 
  ! When setting the jet vs wind flags (jet_wind and jw_switch) use the jet_flag 
  ! and wind_flag parameter values set in the declaration (1 for jets, 2 for winds).
  ! Default to setting flags to 0 if neither
#ifdef debug_jets
  print*, "Particles_wind.F90: Now testing for jets... (w_numloc)", w_numloc
#endif

  ! Test if jets should be on - added 20230726 -SA
  if ( (particles(MASS_PART_PROP, p) .ge. min_jet_mass) .and. &
       (particles(MASS_PART_PROP, p) .lt. max_jet_mass) .and. &
       ( (time - particles(CREATION_TIME_PART_PROP, p)) .lt. jet_time )) then
    ! if in the mass range and less than the age of the jet
    jet_switch = .true.
  
  else 
    ! This could be because there is no feedback or because we want spherical winds
    jet_switch = .false.
  
  end if

  ! Now test for wind condition:  - Added comment and jet_switch 20230726 -SA
  if ( (particles(MASS_PART_PROP, p) .ge. min_wind_mass) .and. & !changed .gt. to .ge. -SA 20230728
     (jet_switch .eqv. .false. ) .and.  &
     (particles(DMDT_PART_PROP, p) .gt. 0.0d0)) then
#ifdef debug_jets
    print*, "Injecting winds for star mass: ", particles(MASS_PART_PROP, p)
#endif
    jet_wind(w_numloc)  = wind_flag  ! Added jet/wind switch -SA 20230726

  ! Now check for jet condition: - Added 20230726 -SA
  else if ( (jet_switch .eqv. .true. ) .and. &
            (particles(DMDT_PART_PROP, p) .gt. 0.0d0)) then
#ifdef debug_jets
    print*, "Injecting jets for star mass: ", particles(MASS_PART_PROP, p)
#endif
    jet_wind(w_numloc)  = jet_flag  ! Added jet/wind switch -SA 20230726

#ifdef debug_jets
  else
    print*, "Particles_wind.F90: Neither jets nor winds turned on..."
#endif 

  end if

  ! print*, "Particles_wind.F90: check jet_wind value: ", jet_wind(w_numloc), jet_flag, wind_flag
  ! print*, "Particles_wind.F90: check particle mass: ", particles(MASS_PART_PROP, p)
  ! call flush()


  if (jet_wind(w_numloc) .gt. 0) then

   ! w_numloc = w_numloc + 1 
  
    locx(w_numloc)      = particles(POSX_PART_PROP, p)
    locy(w_numloc)      = particles(POSY_PART_PROP, p)
    locz(w_numloc)      = particles(POSZ_PART_PROP, p)
    angmom_x(w_numloc)  = particles(X_ANG_PART_PROP, p)
    angmom_y(w_numloc)  = particles(Y_ANG_PART_PROP, p)
    angmom_z(w_numloc)  = particles(Z_ANG_PART_PROP, p) ! Added angular momentum -SA 20230718
    locdmdt(w_numloc)   = particles(DMDT_PART_PROP, p)
    locv_wind(w_numloc) = particles(VELW_PART_PROP, p)
    locc_time(w_numloc) = particles(CREATION_TIME_PART_PROP, p)
    locbgdy(w_numloc)   = particles(BGDY_PART_PROP,p)
    p_ind(w_numloc)     = p 

  end if
end do

! print*, "Particles_wind.F90: test angmom arrays x: ", angmom_x
! print*, "Particles_wind.F90: test angmom arrays y: ", angmom_y
! print*, "Particles_wind.F90: test angmom arrays z: ", angmom_z
! print*, "Particles_wind.F90: test jet wind switch: ", jet_wind

! Now use MPI to vector gather all the information for how to inject
! the winds on each processor.
#ifdef debug
print*, "w_numloc =", w_numloc, dr_globalMe
#endif

disp = 0
rank_minus_one = dr_globalNumProcs - 1

! Gather the array on the root process. Note that we require the
! user to pass the proper length of the final array. This can be
! gotten from get_number_of_new_tags.

! Make an array of the # of incoming particles from each processor.
call MPI_AllGather(w_numloc, 1, MPI_INTEGER, &
	      num_array, 1, MPI_INTEGER, &
	      dr_globalComm, ierr)

! Allocate the actual arrays to pass.

w_num = sum(num_array)

if (allocated(x)) &
    deallocate(x, y, z)
if (allocated(j_x)) &
    deallocate(j_x, j_y, j_z) ! Added angular momentum -SA 20230718
if (allocated(dmdt)) &
    deallocate(dmdt, v_wind, bgdy, c_time)
if (allocated(jw_switch)) &
    deallocate(jw_switch) ! Added 20230726 -SA

allocate(x(w_num), y(w_num), z(w_num))
allocate(j_x(w_num), j_y(w_num), j_z(w_num)) ! Added angular momentum -SA 20230718
allocate(dmdt(w_num), v_wind(w_num), c_time(w_num), bgdy(w_num))
allocate(jw_switch(w_num)) ! Added 20230726 -SA 

x=0.0d0; y=0.0d0; z=0.0d0
j_x=0.0d0; j_y=0.0d0; j_z=0.0d0
dmdt = 0.0d0; v_wind=0.0d0; mass = 0.0d0; c_time = 0.0d0; bgdy=0.0d0
jw_switch = 0

! Set the displacement for the incoming data based on how many
! particles are coming in from each processor. Note the displacement
! for the root process is zero, for rank 1 disp = num on root,
! for rank 2 disp = num on root + num on 1, etc etc.

do i=1, dr_globalNumProcs-1

  disp(i+1) = disp(i) + num_array(i)

end do
#ifdef debug
print*, "About to gather.", dr_globalMe
print*, "num_array =", num_array, dr_globalMe
print*, "disp =", disp, dr_globalMe
#endif
! Now actually gather the info on each proc using the variable length array
! gather command in MPI.
call MPI_AllGatherv(locx, w_numloc, FLASH_REAL, x, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locy, w_numloc, FLASH_REAL, y, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locz, w_numloc, FLASH_REAL, z, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(angmom_x, w_numloc, FLASH_REAL, j_x, num_array, &
           disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(angmom_y, w_numloc, FLASH_REAL, j_y, num_array, &
           disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(angmom_z, w_numloc, FLASH_REAL, j_z, num_array, &
           disp, FLASH_REAL, dr_globalComm, ierr)  ! Added angular momentum - SA 20230718
call MPI_AllGatherv(locdmdt, w_numloc, FLASH_REAL, dmdt, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locv_wind, w_numloc, FLASH_REAL, v_wind, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locc_time, w_numloc, FLASH_REAL, c_time, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locbgdy, w_numloc, FLASH_REAL, bgdy, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(jet_wind, w_numloc, MPI_INTEGER, jw_switch, num_array, &
           disp, MPI_INTEGER, dr_globalComm, ierr)  ! Fixed var type - SA 20230728

!!!!! TODO: Consider combining all of these into a single MPI_AllGatherv with a large 2D array containing all properties

! Now all procs have an array of each value in the same order, so we can
! inject the wind at each point across all procs.
#ifdef debug
print*, "Done gathering.", dr_globalMe
#endif

! print*, "Starting loop with inject_direct call (w_num): ", w_num , " -SA 202212"
! print*, "Before inject_direct - ang_mom:", [j_x, j_y, j_z], [angmom_x, angmom_y, angmom_z]

do p=1, w_num
  !dmdt(p) = 1d-6*solarMass/yr
  mass  = dmdt(p)*dt ! Total mass injected by this star this step.
  twind = dr_simTime + dt - c_time(p) ! Time since the start of this stars wind.
  bgdy_old = bgdy(p) ! Background density of the gas when the wind started.
#ifdef debug2
  if (dr_globalMe .eq. 0) then
    print*, "Calling inject direct with mass, dt, dmdt, vwind, bgdy =", mass, dt, dmdt(p)/solarMass*yr, v_wind(p), bgdy(p)
    print*, "Calling inject direct with angular momentum vector and jet/wind: ", [j_x(p), j_y(p), j_z(p)], jw_switch(p) !-SA 202307
    print*, "index of loop: ", p, "and position: ", x(p), y(p), z(p), " -SA 202212"
  endif
#endif
  
  call inject_direct([x(p), y(p), z(p)], [j_x(p), j_y(p), j_z(p)], jw_switch(p), mass, v_wind(p), mass, twind, dt, bgdy(p)) 
  !Added j_i -SA 20230718
  !Added jw_switch -SA 20230726  Fix order 20230802 SA

! If this call to inject_direct calculated the background density, store it on the proper processor.
  if (bgdy_old .eq. 0.0d0) then ! no recorded background density, so must be first loop.
    if (disp(dr_globalMe+1) - p < 0) then ! do I own this particle?
      if (dr_globalMe .eq. dr_globalNumProcs - 1) then
        particles(BGDY_PART_PROP, p_ind(p-disp(dr_globalMe+1))) = bgdy(p)
      else if (disp(dr_globalMe+2) - p >= 0) then
        particles(BGDY_PART_PROP, p_ind(p-disp(dr_globalMe+1))) = bgdy(p)
      end if
    end if
  end if

end do

! print*, "Particles_wind.F90: check deallocations: ", jet_wind, jw_switch

! Add a flush call before deallocation so we can check what's going on
! print*, "Particles_wind.F90: Implementing flush call now."
! call flush()

deallocate(p_ind)

deallocate(locx, locy, locz)
deallocate(angmom_x, angmom_y, angmom_z) ! Add ang momentum -SA 20230720
deallocate(locdmdt, locv_wind, locbgdy, locc_time)
deallocate(dmdt, v_wind, c_time, bgdy)
deallocate(x, y, z)
deallocate(j_x, j_y, j_z)  ! Added ang momentum -SA 20230720 
deallocate(jet_wind, jw_switch) !Added -SA 20230726
! Let the Grid unit know we updated these variables to properly fill guard cells.

call Grid_notifySolnDataUpdate() !(/ EINT_VAR, ENER_VAR, TEMP_VAR, VELX_VAR, VELY_VAR, VELZ_VAR, DENS_VAR /)

call Grid_fillGuardCells(CENTER, ALLDIR) !, doEos=.true., eosMode=MODE_DENS_EI, selectBlockType=ACTIVE_BLKS)

end subroutine Particles_wind

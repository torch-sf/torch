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

!#define debug
#define debug2


subroutine Particles_wind(dt)

use Particles_data, only : particles, pt_typeInfo, pt_numLocal
use Particles_windData
use Particles_interface, only : Particles_getGlobalNum

!use Particles_windData

use Driver_interface, only : Driver_abortFlash
use Driver_data, only : dr_globalComm, dr_globalNumProcs, dr_globalMe, dr_simTime

use Grid_interface, only : Grid_fillGuardCells

use RuntimeParameters_interface, ONLY: RuntimeParameters_get

! Timers used to time entire Particle_winds, MPI AllGather calls and inject direct routine.
use Timers_interface, ONLY : Timers_start, Timers_stop  !SA 20240207


implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"
#include "Particles.h"
#include "Eos.h"

real, intent(in)       :: dt

! Injection info for winds.
! Global vars
real*8, allocatable :: x(:), y(:), z(:), &
               dmdt(:), v_wind(:), c_time(:), bgdy(:)

integer    :: w_num
! Local vars
real, allocatable      :: locx(:), locy(:), locz(:), locdmdt(:), locv_wind(:), &
			  locc_time(:), locbgdy(:)

real                   :: mass, twind, bgdy_old

! For counting particles.
integer                :: p_begin, p_end, p_num, p_globalnum, w_numloc

! Storage of the local indices that are actual wind stars.
integer, allocatable   :: p_ind(:)

! For MPI comm
integer                :: num_array(dr_globalNumProcs), &
                          disp(dr_globalNumProcs), ierr, i, p, &
			  rank_minus_one

! Runtime parameters
character(len=8), save :: part_type
real, parameter :: yr = (60.0d0**2.0)*24.0d0*365.25d0
real, parameter :: solarMass = 1.989d33
logical, save          :: first_call = .true.

#ifdef TRACER_FIELDS
! For tracer fields
real*8, allocatable           :: dydt(:,:)
real, allocatable             :: locdydt(:,:)
integer, parameter            :: pt_num_tracer_fields = NMASS_SCALARS
integer, parameter            :: pt_tracer_fields_begin = Y001_PART_PROP
real, dimension(pt_num_tracer_fields) :: yields
integer, dimension(pt_num_tracer_fields), save :: pt_tracer_inds
integer :: itracer
#endif

if (first_call) then
  call RuntimeParameters_get("min_wind_mass", min_wind_mass)
  
#ifdef TRACER_FIELDS
  ! Gather tracer field indices for particle yields.
  do itracer = 1, pt_num_tracer_fields
    pt_tracer_inds(itracer) = pt_tracer_fields_begin + (itracer - 1)
  enddo
#endif

  first_call = .false.
end if

min_wind_dt = 1d99

call Timers_start("Particles_wind")

! Local number of massive/active particles.
p_begin = pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE)
p_num   = pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE)
p_end   = p_num + p_begin - 1

allocate(p_ind(pt_numLocal))
p_ind = 0

call Particles_getGlobalNum(p_globalnum)

! First, we allocate the length of each array to be the total number of star particles.
! Each array is then initialized to be all 0. -SA 20240216

allocate(locx(p_globalnum), locy(p_globalnum), locz(p_globalnum), locc_time(p_globalnum))
allocate(locdmdt(p_globalnum), locv_wind(p_globalnum), locbgdy(p_globalnum))

w_numloc  = 0
w_num     = 0
num_array = 0

locx = 0.0d0; locy=0.0d0; locz=0.0d0
locdmdt = 0.0d0; locv_wind=0.0d0; locbgdy=0.0d0; locc_time= 0.0d0

#ifdef TRACER_FIELDS
allocate(locdydt(p_globalnum,pt_num_tracer_fields))
locdydt = 0.0d0
#endif

! This do loop then loops over each particle and checks if winds are on for that star.
! If winds are on, then a separate index w_numloc (which starts at 0) is incremented and
! the entry of each of the above arrays which corresponds to the new w_numloc index is
! set to the corresponding non-zero value for the star. At the end of the loop the first
! entries (a number matching the number of wind stars) will be non-zero and all subsequent
! entries will be zero.  The final value of the w_numloc index will also track the
! total number of wind stars. -SA 20240216

do p = p_begin, p_end
#ifdef debug
  print*, "Particle mass =", particles(MASS_PART_PROP, p)
  print*, "Particle dmdt =", particles(DMDT_PART_PROP, p)
#endif
  if ( (particles(MASS_PART_PROP, p) .gt. min_wind_mass) .and. &
     (particles(DMDT_PART_PROP, p) .gt. 0.0d0)) then
     
    w_numloc = w_numloc + 1
  
    locx(w_numloc)      = particles(POSX_PART_PROP, p)
    locy(w_numloc)      = particles(POSY_PART_PROP, p)
    locz(w_numloc)      = particles(POSZ_PART_PROP, p)
    locdmdt(w_numloc)   = particles(DMDT_PART_PROP, p)
    locv_wind(w_numloc) = particles(VELW_PART_PROP, p)
    locc_time(w_numloc) = particles(CREATION_TIME_PART_PROP, p)
    locbgdy(w_numloc)   = particles(BGDY_PART_PROP,p)
    p_ind(w_numloc)     = p

#ifdef TRACER_FIELDS
    do itracer = 1,pt_num_tracer_fields
      locdydt(w_numloc,itracer) = particles(pt_tracer_inds(itracer), p)
    enddo
#endif

  end if
end do

! Now that the above arrays have identified all the wind stars for a
! given processor, the following code collects that info from all
! processors using MPI_AllGather. Only the first w_numloc entries
! from each processor are gathered, meaning all the zero valued
! entries should be dropped during the MPI_AllGather stage. The
! final arrays should only have wind stars with non-zero values. -SA 20240216

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
if (allocated(dmdt)) &
    deallocate(dmdt, v_wind, bgdy, c_time)

allocate(x(w_num), y(w_num), z(w_num))
allocate(dmdt(w_num), v_wind(w_num), c_time(w_num), bgdy(w_num))

x=0.0d0; y=0.0d0; z=0.0d0
dmdt = 0.0d0; v_wind=0.0d0; mass = 0.0d0; c_time = 0.0d0; bgdy=0.0d0

#ifdef TRACER_FIELDS
! Reallocate and reset arrays (safety)
if (allocated(dydt)) deallocate(dydt)
allocate(dydt(w_num,pt_num_tracer_fields))
yields = 0.0d0
dydt = 0.0d0
#endif

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

call Timers_start("MPI_AllGather_winds")

! Now actually gather the info on each proc using the variable length array
! gather command in MPI.
call MPI_AllGatherv(locx, w_numloc, FLASH_REAL, x, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locy, w_numloc, FLASH_REAL, y, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locz, w_numloc, FLASH_REAL, z, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locdmdt, w_numloc, FLASH_REAL, dmdt, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locv_wind, w_numloc, FLASH_REAL, v_wind, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locc_time, w_numloc, FLASH_REAL, c_time, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locbgdy, w_numloc, FLASH_REAL, bgdy, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)

#ifdef TRACER_FIELDS
do itracer = 1,pt_num_tracer_fields
   call MPI_AllGatherv(locdydt(:,itracer), w_numloc, FLASH_REAL, dydt(:,itracer), num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
enddo
#endif

call Timers_stop("MPI_AllGather_winds")

! Now all procs have an array of each value in the same order, so we can
! inject the wind at each point across all procs.
#ifdef debug
print*, "Done gathering.", dr_globalMe
#endif

! The following do loop goes over all the entries of the gathered arrays and injects
! winds with inject-direct. At this point every entry of these arrays should be
! a wind star with non-zero dmdt. Thus, w_num is the number of wind stars. -SA 20240216

do p=1, w_num
  !dmdt(p) = 1d-6*solarMass/yr
  mass  = dmdt(p)*dt ! Total mass injected by this star this step.
  twind = dr_simTime + dt - c_time(p) ! Time since the start of this stars wind.
  bgdy_old = bgdy(p) ! Background density of the gas when the wind started.

#ifdef TRACER_FIELDS
  do itracer = 1,pt_num_tracer_fields
    yields(itracer) = dydt(p,itracer)*dt ! Total mass of this tracer field injected by this star this step.
  enddo
#endif

#ifdef debug2
    if (dr_globalMe .eq. 0) &
      print*, "Calling inject direct with inj mass, dt, dmdt, vwind, bgdy =", mass, dt, dmdt(p)/solarMass*yr, v_wind(p), bgdy(p)
#endif

  call Timers_start("inject_direct_call")
#ifdef TRACER_FIELDS
  call inject_direct([x(p), y(p), z(p)], mass, v_wind(p), yields, twind, dt, bgdy(p))
#else
  call inject_direct([x(p), y(p), z(p)], mass, v_wind(p), 0.0, twind, dt, bgdy(p)) !Remove duplicate mass -SA 20240207
#endif
  call Timers_stop("inject_direct_call")

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

deallocate(p_ind)

deallocate(locx, locy, locz)
deallocate(locdmdt, locv_wind, locbgdy, locc_time)
deallocate(dmdt, v_wind, c_time, bgdy)
deallocate(x, y, z)

#ifdef TRACER_FIELDS
deallocate(locdydt)
deallocate(dydt)
#endif

! Let the Grid unit know we updated these variables to properly fill guard cells.

call Grid_notifySolnDataUpdate() !(/ EINT_VAR, ENER_VAR, TEMP_VAR, VELX_VAR, VELY_VAR, VELZ_VAR, DENS_VAR /)

call Grid_fillGuardCells(CENTER, ALLDIR) !, doEos=.true., eosMode=MODE_DENS_EI, selectBlockType=ACTIVE_BLKS)

call Timers_stop("Particles_wind")

end subroutine Particles_wind

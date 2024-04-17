!!****f* source/physics/RadTrans/RadTrans_computeDt
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  RadTrans_computeDt(integer(IN) :: blockID,
!!                     integer(IN) :: blkLimits(2,MDIM)
!!                     integer(IN) :: blkLimitsGC(2,MDIM)
!!                     real(IN),pointer::  solnData(:,:,:,:),   
!!                     real(OUT)   :: dt_radtrans, 
!!                     real(OUT)   :: dt_minloc(5)) 
!!  DESCRIPTION 
!!    Compute radiative transfer time step
!!
!!  ARGUMENTS
!!    blockID       --  local block ID
!!    blkLimits     --  the indices for the interior endpoints of the block
!!    blkLimitsGC   --  the indices for endpoints including the guardcells
!!    solnData      --  the physical, solution data from grid
!!    dt_radtrans   --  variable to hold timestep constraint
!!    dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!
!!***
! DORIC does subcycling so no need for this


!#define DEBUG_RADDT
!#define ONLY_FB_PROC

subroutine RadTrans_computeDt(blockID,  blkLimits,blkLimitsGC, &
     solnData, dt_radtrans, dt_minloc)

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  use Driver_data, only : dr_globalMe, dr_nstep, dr_globalComm, dr_simTime

#ifdef WIND_INJ
  use Particles_windData, only : min_wind_mass, min_wind_dt, &
#ifdef JETS
                                 min_jet_mass, max_jet_mass, jet_time, & ! -SA 20230920
#endif
                                 wind_target_temp

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
#endif

  use Particles_rayData, only : ph_radPressure, cfl_radPressure
  use rt_data, only: rt_dt, rt_dt_pos, rt_protonMass, rt_gamma1, &
                     rt_rayTrace, rt_useNumstepsRadTransDtOnStart, &
                     rt_numstepsRadTransDt, rt_useRadTransDt, rt_dt_temp
  use Particles_data, only: particles, pt_typeInfo
  use Grid_interface, only: Grid_getMinCellSize
  use Particles_windData, ONLY: mass_load, wind_target_temp
  implicit none
  
#include "Flash_mpi.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: blkLimits(2,MDIM)
  integer, intent(IN) :: blkLimitsGC(2,MDIM)
  real, pointer :: solnData(:,:,:,:) 
  real, intent(INOUT) :: dt_radtrans
  integer, intent(INOUT)  :: dt_minloc(5)
  integer, save       :: num_fb_stars=0 ! number of feedback producing stars.
  integer             :: p, mass_count, jet_mass_count, type_begin, type_count, type_end, ierr
  real, parameter     :: eightMSun = 8.0d0*1.989d33, kB = 1.381d-16, mu = 0.61
  real                :: deltaX, csRad, photon_flux, max_mass, max_mass_jet, ener_per_ph, cross_sec
  real                :: dt_mom, dt_wind, dt_jet, wind_vel, jet_vel, dt_min_local
  ! This makes the code use the RadTransDt for a set number of loops after
  ! a new feedback star appears.
  !integer, parameter  :: maxnum_radTransDt_loops = 100
  integer, save       :: currnum_radTransDt_loops = 999, curr_step_num = -1
  integer, save       :: global_mass_count=0, global_fb_stars=0
  logical, save       :: stillUseRadDt = .false., first_call=.true.
  real, save          :: old_dt_radtrans = 1d99
  real :: mass_load_factor
  real :: injectVelocity
  real :: time ! Current simulation time for determining star particle age -SA 20240226
  character(len=10), save :: conserved_quant
  real :: refVel
  real, save :: sink_density

  ! Here we estimate the timestep for feedback in radiation and winds
  ! from a star we have introduced to the simulation in between loop steps.
  ! Note this is necessary (and also to have the global dt calculated 
  ! before any call to Hydro) to ensure that the Hydro knows about the
  ! proper velocities and sound speed and doesn't take too large a step. - JW

  ! Also since this makes worst case assumptions, we don't actually use any
  ! grid information. So this only needs to run once per processor. Therefore
  ! after one run, we just return on all other calls. - JW

#ifdef WIND_INJ
  if (first_call) then
    call RuntimeParameters_get("cons_quant", conserved_quant)
    call RuntimeParameters_get("min_wind_mass", min_wind_mass)
#ifdef JETS
    call RuntimeParameters_get("min_jet_mass", min_jet_mass) ! -SA 20230920
    call RuntimeParameters_get("max_jet_mass", max_jet_mass)
    call RuntimeParameters_get("jet_time", jet_time)
#endif
    call RuntimeParameters_get("wind_target_temp", wind_target_temp)
    call RuntimeParameters_get("rt_useRadTransDt", rt_useRadTransDt)
    call RuntimeParameters_get("rt_useNumstepsRadTransDtOnStart", rt_useNumstepsRadTransDtOnStart)
    call RuntimeParameters_get("rt_numstepsRadTransDt", rt_numstepsRadTransDt)
    call RuntimeParameters_get("rt_dt_temp", rt_dt_temp)
    call RuntimeParameters_get("cfl_radPressure", cfl_radPressure)
    call RuntimeParameters_get("sink_density_thresh", sink_density)
    first_call = .false.
  end if

#ifdef JETS
  !!! Note JETS options assume WIND_INJ also is on -SA 20240223
  ! If checking for jets, we need the current simulation time -SA 20240223
  time = dr_simTime !no modification by dt
#endif
#endif

  if (.not. rt_useRadTransDt) then
    dt_radtrans = 1e99
    return
  end if


  mass_count = 0
  jet_mass_count = 0
  max_mass   = 0.0d0
  max_mass_jet = 0.0d0
  csRad      = 0.0d0
  dt_mom     = 1d99
  dt_wind    = 1d99
  dt_jet     = 1d99
  wind_vel   = 0.0d0
  jet_vel    = 0.0d0
  dt_min_local = 1d99

  if (curr_step_num /= dr_nstep) then ! We have a new dt calc loop starting
    curr_step_num = dr_nstep ! initialize current step
  else
    return ! Otherwise immediately return.
  end if

  type_begin = pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE)
  type_count = pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE)
  type_end   = type_count + type_begin - 1

  call Grid_getMinCellSize(deltaX)
  
#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
  if (stillUseRadDt) then
#endif
  print*, "On entry.", dr_globalMe
  print*, "currnum_radTransDt_loops=", currnum_radTransDt_loops, dr_globalMe
  print*, "stillUseRadDt=", stillUseRadDt, dr_globalMe
  print*, "min_wind_mass=", min_wind_mass, dr_globalMe
#ifdef JETS
  print*, "min_jet_mass=",  min_jet_mass, dr_globalMe ! -SA 20230920
#endif
#ifdef ONLY_FB_PROC
  endif
#endif
#endif
  
  do p=type_begin, type_end
#if defined(JETS)
    ! If JETS is on, we need to make sure we have a wind and NOT a jet star
    if (particles(MASS_PART_PROP, p) .ge. min_wind_mass) then !if producing winds
       if ((particles(MASS_PART_PROP, p) .lt. min_jet_mass) .or. & !less massive than jet range
           (particles(MASS_PART_PROP, p) .ge. max_jet_mass) .or. & !more massive than jet range
           ( (time - particles(CREATION_TIME_PART_PROP, p)) .ge. jet_time )) then !older than jets
#elif defined(WIND_INJ)
    ! For wind dt, just use most massive star and wind parameters -SA 20240223
    if (particles(MASS_PART_PROP, p) .ge. min_wind_mass) then
#else
    ! TODO: We should eventually update this to not hardcode a mass value -SA 20230920
    if (particles(MASS_PART_PROP, p) .ge. eightMSun) then
#endif
      mass_count = mass_count + 1
      if (particles(MASS_PART_PROP, p) .gt. max_mass) then ! We want the most massive star flux.
        max_mass    = particles(MASS_PART_PROP, p)
        photon_flux = particles(NION_PART_PROP, p) ! # of photons per sec from the star.
        ener_per_ph = particles(EION_PART_PROP, p) ! Average energy per photon.
        cross_sec   = particles(SIGH_PART_PROP, p) ! Cross section of hydrogen.
#ifdef WIND_INJ
        wind_vel    = particles(VELW_PART_PROP, p) ! Wind velocity.
#endif
        photon_flux = photon_flux ! # Assume worst case, star dumps all into one cell.
      end if
#ifdef JETS
      end if
#endif
    end if

#ifdef JETS
    !!  Now look for most massive jet star to calculate det_jet -SA 20240223
    if ((particles(MASS_PART_PROP, p) .ge. min_jet_mass) .and. &
        (particles(MASS_PART_PROP, p) .lt. max_jet_mass) .and. &
        ( (time - particles(CREATION_TIME_PART_PROP, p)) .lt. jet_time )) then
       jet_mass_count = jet_mass_count + 1
       if (particles(MASS_PART_PROP, p) .gt. max_mass_jet) then ! We want the most massive jet star actually spawned.
           max_mass_jet    = particles(MASS_PART_PROP, p)
           jet_vel    = particles(VELW_PART_PROP, p) ! Jet velocity - same param. as wind_vel but for jet star
       end if
    end if
#endif
  end do

#ifdef JETS
  !! Note: This will double count any stars that currently are producing jets but
  !! will eventually produce winds. There might be a better way to do this. -SA 20240223
  mass_count = mass_count + jet_mass_count
#endif

  ! Get the total mass count.
  call MPI_ALLREDUCE(mass_count, global_mass_count, 1, MPI_INT, MPI_SUM, &
                     dr_globalComm, ierr)

#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
  if (stillUseRadDt) then
#endif
    print*, "Before check.", dr_globalMe
    print*, "currnum_radTransDt_loops=", currnum_radTransDt_loops, dr_globalMe
    print*, "stillUseRadDt=", stillUseRadDt, dr_globalMe
    print*, "mass_count, num_fb_stars=", mass_count, num_fb_stars, dr_globalMe
    print*, "global_mass_count, global_fb_stars=", global_mass_count, global_fb_stars, dr_globalMe
    print*, "max_mass=", max_mass/1.989e33, dr_globalMe
#ifdef ONLY_FB_PROC
  end if
#endif
#endif
  ! If both 1) We have a new FB star that hasn't been on the grid before and
  !         2) Your the proc that currently has it.
  ! Then compute an estimate of the timestep.
  !! I am concerned that since mass_count can go down if stars leave the grid (or a jet reaches the end
  !! its life) this check isn't very meaningful. TODO is there a better way to do this? -SA 20240223
  if ((global_mass_count .gt. global_fb_stars) .and. (mass_count .gt. num_fb_stars)) then ! We have a new feedback star.
    currnum_radTransDt_loops = 0
  end if

! If we haven't gone maxnum steps, we should still be using the RadTransDt.
  stillUseRadDt =  (currnum_radTransDt_loops .lt. rt_numstepsRadTransDt)

  if (stillUseRadDt) then 
    currnum_radTransDt_loops = currnum_radTransDt_loops + 1

#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
    if (stillUseRadDt) then
#endif
      print*, "After check.", dr_globalMe
      print*, "currnum_radTransDt_loops=", currnum_radTransDt_loops, dr_globalMe
      print*, "stillUseRadDt=", stillUseRadDt, dr_globalMe
      print*, "mass_count, num_fb_stars=", mass_count, num_fb_stars, dr_globalMe
      print*, "max_mass=", max_mass/1.989e33, dr_globalMe
#ifdef ONLY_FB_PROC
    end if
#endif
#endif
    
    if (currnum_radTransDt_loops == 1) then
#ifdef WIND_INJ
      if (wind_vel > 0.0) then
! Assume new star winds will make the gas ionized and 2x10^7 K. -JW
! This temperature should be set to the wind_target_temp of mass loading. -BP 06.23.22
        csRad        = sqrt(rt_gamma1 * kB * wind_target_temp / rt_protonMass / mu)

! The velocity used for wind timestep should account for mass loading!!
if (mass_load) then
    refVel = sqrt(wind_target_temp/1.38d7)*1e8
    if (conserved_quant .eq. "momentum") then
        mass_load_factor = wind_vel/refVel - 1.0d0
        injectVelocity   = wind_vel / (1.0d0+mass_load_factor)
    else if (conserved_quant .eq. "energy") then
        mass_load_factor = wind_vel**2.0d0/refVel**2.0d0 - 1.0d0
        injectVelocity   = wind_vel / dsqrt(1.0d0+mass_load_factor)
    endif
    ! Move this within if statement -SA 20240223
    dt_wind      = 0.3 * deltaX / sqrt(csRad**2.0 + injectVelocity**2.0)

end if

#ifdef JETS
        !!! Add new jets dt - just uses Courant condition. The 0.3 factor is empirical
        !!! (matching dt_wind, above) and could be improved on in the future. -SA 20240417
        dt_jet = 0.3 * deltaX / jet_vel
#endif

#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
        if (stillUseRadDt) then
#endif
          write(*,'(A,X,ES13.3,I4)') "csRad after wind =", csRad, dr_globalMe
          write(*,'(A,X,ES13.3,I4)') "dt_wind =", dt_wind, dr_globalMe
#ifdef ONLY_FB_PROC
        end if
#endif
#endif
      end if
#endif

      if (rt_rayTrace) then
! Assume new star radiation will make the gas ionized and 2x10^5 K. -JW
! This should be a flash runtime parameter (rt_dt_temp). -BP 06.23.22
        csRad        = sqrt(rt_gamma1 * kB * rt_dt_temp / rt_protonMass / mu)
    
! Estimate the velocity from the momentum of radiation as that which
! makes it into one of the surrounding cells based on the newly formed
! stars photon flux on that cell. Note this is a bit too aggressive, should
! go back and actually integrate and find the flux to properly put in here. - JW

! It's ok to assume that in the initial timestep all the ionizing photons get
! absorbed because the sink density in which a star will form implies a mean 
! free path smaller than a cell width. -BP 19Jan23

        if (ph_radPressure) then
! This is delta_t = sqrt(V * c * mu * m_H / (dN_ph/dt * sigma_H * E_avg) - JW
!          dt_mom        = 0.3*sqrt(deltaX**3.0d0 * 3.0d10 * 1.24d0 * 1.67d-24 / &
!                        (photon_flux * cross_sec * ener_per_ph))

! Above formula is dt = mean free path / dv, when it should be dt = dx/dv 
! where dv = (E_avg * dN_ph/dt * dt) / (rho * V * c). Density is approximated as
! the sink threshhold density, as this is roughly the density where newborn stars 
! will form in. Also made the "courant" factor a user parameter. -BP 19Jan23
           dt_mom = cfl_radPressure * sqrt(deltaX**4.0d0 * 3.0d10 * sink_density / &
                               (photon_flux * ener_per_ph))

        end if
      dt_min_local  = 0.3*deltaX / csRad ! With a safe CFL factor for now.

#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
        if (stillUseRadDt) then
#endif
          write(*,'(A,X,ES13.3,I4)') "csRad after rad =", csRad, dr_globalMe
          write(*,'(A,X,ES13.3,I4)') "dt_mom =", dt_mom, dr_globalMe
          write(*,'(A,X,ES13.3,I4)') "dt_radtrans =", dt_radtrans, dr_globalMe
#ifdef ONLY_FB_PROC
        end if
#endif
#endif    
      end if
!!! Add jets dt to determination of min dt value -SA 20240223
      dt_min_local  = min(dt_min_local, dt_wind, dt_mom, dt_jet)
      old_dt_radtrans = dt_min_local
#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
      if (stillUseRadDt) then
#endif
        print*, "RadTrans dt = ", dt_min_local, dr_globalMe
#ifdef ONLY_FB_PROC
      end if
#endif
#endif
    else
      dt_min_local = old_dt_radtrans
    end if
  else ! Hydro should already have the proper dt any other time.
    currnum_radTransDt_loops = 999
#ifdef WIND_INJ
    dt_min_local = min_wind_dt
#else
    dt_min_local  = 1e99 !rt_dt
#endif
  end if
  
  if (dt_min_local .lt. dt_radtrans) dt_radtrans = dt_min_local
#ifdef DEBUG_RADDT
#ifdef ONLY_FB_PROC
    if (stillUseRadDt) then
#endif
      write(*,'(A,X,ES13.3,I4)') "dt_radtrans at end=", dt_radtrans, dr_globalMe
#ifdef ONLY_FB_PROC
    end if
#endif
#endif 
  
  if (dt_radtrans .le. 0.0) dt_radtrans = 1d99

  dt_minloc      = rt_dt_pos

  ! Update the total numbers of feedback stars. Note this number can
  ! go down if a massive star leaves the system.
  num_fb_stars = mass_count
  call MPI_ALLREDUCE(num_fb_stars, global_fb_stars, 1, MPI_INT, MPI_SUM, &
                     dr_globalComm, ierr)
  
  return
end subroutine RadTrans_computeDt

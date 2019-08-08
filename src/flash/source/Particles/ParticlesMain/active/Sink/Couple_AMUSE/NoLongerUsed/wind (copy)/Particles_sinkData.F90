!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkData
!!
!! NAME
!!
!!    Particles_sinkData
!!
!! SYNOPSIS
!!
!!    Particles_sinkData()
!!
!! DESCRIPTION
!!
!!    Module to hold local variables and data types for sink particle unit
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!***

module Particles_sinkData

  implicit none

#include "constants.h"
#include "Particles.h"
#include "Flash.h"

  public :: MAX_MSGS, maxsinks, n_empty
  logical, save :: RunningParticles = .false.
  integer, save :: n_empty
  integer, save :: sink_maxSinks
  logical, save :: useSinkParticles
  logical, save :: sink_AdvanceSerialComputation

  integer, parameter :: maxsinks = 20000 !1024
  integer, save, pointer, dimension(:)  :: NumParticlesPerBlock
  integer, parameter :: MAX_MSGS = 12

  real, dimension(MAX_MSGS) :: send_buff, recv_buff

  integer, save :: ipx, ipy, ipz, ipvx, ipvy, ipvz, ipm, iptag
  integer, save :: ipblk, iplx, iply, iplz, iplx_old, iply_old, iplz_old, ipmdot, ipt
  integer, save :: ipdtold, ipcpu, iold_pmass

  integer, parameter :: pt_sinkParticleProps = NPART_PROPS

  ! particles_local and particles_global refer to 
  ! sink particles - the local list and the global list
  real, save, allocatable, dimension(:,:) :: particles_local
  real, save, allocatable, dimension(:,:) :: particles_global
  
  integer, save :: local_tag_number
  integer, save :: localnp, localnpf

  ! variables needed for Ewald correction in case of periodic BCs
  integer, save :: sink_EwaldNx, sink_EwaldNy, sink_EwaldNz
  integer, save :: sink_EwaldSeriesN
  character(len=MAX_STRING_LENGTH), save :: sink_EwaldFileName
  real, save :: sink_EwaldDxI, sink_EwaldDyI, sink_EwaldDzI
  real, save, allocatable, dimension(:,:,:) :: sink_EwaldFieldX, &
                                               sink_EwaldFieldY, &
                                               sink_EwaldFieldZ

  ! Variables for AMUSE to learn if new sinks were formed during a
  ! Flash evolution step. - JW
  
  integer, save :: number_new_sinks
  integer, save, allocatable, dimension(:) :: new_tags
  logical, save :: new_evolution_step = .true.
  real, save :: min_wind_dt = 1e99

end module

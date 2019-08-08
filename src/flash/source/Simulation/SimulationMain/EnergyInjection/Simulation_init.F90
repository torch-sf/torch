!!****if* source/Simulation/SimulationMain/EnergyInjection/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!  Initializes initial conditions for EnergyInjection test problem
!!
!! ARGUMENTS
!!
!!   
!!
!!
!!***
subroutine Simulation_init() 
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use pt_sourceUtil
  
  use Driver_data, ONLY : dr_globalMe
  use Particles_sinkData, ONLY : particles_local
  use pt_sinkInterface, ONLY : pt_sinkCreateParticle, pt_sinkGatherGlobal
  use rt_data, ONLY : eth0, ev2erg, ah0, rt_vary_atomic_frac
  use Particles_interface, ONLY : Particles_sinkMoveParticles, &
      Particles_sinkSortParticles
  implicit none

  integer :: blockID, pno ! blockID, particle #
  real    :: pt           ! particle creation time
  logical :: manual_sink_add

#include "Flash.h"
#include "constants.h"

  ! get the runtime parameters relevant for this problem

  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallX', sim_smallX)

  call RuntimeParameters_get('particle_1_x', sim_p1x)
  call RuntimeParameters_get('particle_1_y', sim_p1y)
  call RuntimeParameters_get('particle_1_z', sim_p1z)

  call RuntimeParameters_get('num_sources', sim_nPtot)

  call RuntimeParameters_get('eos_singleSpeciesA', sim_molarMass)
  call PhysicalConstants_get('ideal gas constant', sim_gasconstant)
  call PhysicalConstants_get('proton mass', sim_protonMass)

  call RuntimeParameters_get('amTemp', sim_amTemp)
  call RuntimeParameters_get('amNumDens', sim_amNumDens)

  call RuntimeParameters_get('sim_Eph', sim_Eph)
  call RuntimeParameters_get('sim_Nph', sim_Nph)
  call RuntimeParameters_get('sim_init_Hp', sim_init_Hp)

#ifdef VARY_ATM_FRAC
  call RuntimeParameters_get('rt_vary_atomic_frac', rt_vary_atomic_frac)
#endif
  
  call RuntimeParameters_get( 'bx0', sim_magx)
  call RuntimeParameters_get( 'by0', sim_magy)
  call RuntimeParameters_get( 'bz0', sim_magz)

  call RuntimeParameters_get( 'sim_tdust', sim_tdust)

  sim_abar = 1.0 + sim_abundM*sim_metal
  
! Now I'll drop in a sink particle. -JW

manual_sink_add = .false.
if (manual_sink_add) then
    if (dr_globalMe == MASTER_PE) then

        blockID = 1
        pt = 0.0

        pno = pt_sinkCreateParticle(sim_p1x, sim_p1y, sim_p1z, pt, blockID, dr_globalMe)

        particles_local(VELX_PART_PROP, 1) = 0.0
        particles_local(VELY_PART_PROP, 1) = 0.0
        particles_local(VELZ_PART_PROP, 1) = 0.0
        particles_local(MASS_PART_PROP, 1) = 30.0*1.989e33 ! 8 solar masses, for now.
        particles_local(NION_PART_PROP, 1) = sim_Nph
        particles_local(EION_PART_PROP, 1) = sim_Eph*ev2erg
        particles_local(SIGH_PART_PROP, 1) = 2.6e-18 !ah0

! Testing for photoelectric heating. - JW
#ifdef PE_HEAT
        particles_local(NPEP_PART_PROP, 1) = 10**(48.6)
        particles_local(EPEP_PART_PROP, 1) = 8.78*ev2erg
        particles_local(SPEP_PART_PROP, 1) = 1.34d-21
#endif

        write(*,'(A,4(1X,ES16.9),3I8)') "initial sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
          & sim_p1x, sim_p1y, sim_p1z, pt, blockID, sim_meshMe, int(particles_local(TAG_PART_PROP,pno))

    endif

    call pt_sinkGatherGlobal()
end if

end

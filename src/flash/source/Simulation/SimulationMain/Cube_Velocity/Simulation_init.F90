!!****if* source/Simulation/SimulationMain/Cube_Velocity/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer myPE)
!!
!! ARGUMENTS
!!
!!    myPE      Current Processor Number
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get routine for initialization.
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getComm

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: i,j,k, istat, ii,jj,kk
  character(len=255),save :: dumstr

  real :: eint, ek ! for sim_cubeFile velocities
  integer :: sim_comm, sim_myPE


  !AT 20190221 - use the multispecies gamma instead.  We are always using 5/3 so no effect
  !call RuntimeParameters_get('gamma', sim_gamma) ! overlap w/ Eos_init.F90
  call RuntimeParameters_get('smlrho', sim_smlrho) ! declared by Hydro/HydroMain/unsplit
  call RuntimeParameters_get('smallp', sim_smallp) ! declared by Hydro/HydroMain/unsplit
  !call RuntimeParameters_get('smallX', sim_smallX) ! declared by Hydro/HydroMain/unsplit

  call RuntimeParameters_get('sim_bx0', sim_magx)
  call RuntimeParameters_get('sim_by0', sim_magy)
  call RuntimeParameters_get('sim_bz0', sim_magz)

  ! stratbox initialization
  call RuntimeParameters_get('sim_useStrat',  sim_useStrat)
  call RuntimeParameters_get('sim_p',  sim_p)
  call RuntimeParameters_get('sim_rho',  sim_rho)
  call RuntimeParameters_get('sim_pIGM',  sim_pIGM)
  call RuntimeParameters_get('sim_rhoIGM',  sim_rhoIGM)
  ! stratbox, gravity profile parameters, stellar disk
  call RuntimeParameters_get('sim_aParm1',  sim_aParm1)
  call RuntimeParameters_get('sim_aParm2',  sim_aParm2)
  call RuntimeParameters_get('sim_aParm3',  sim_aParm3)
  call RuntimeParameters_get('sim_aParm4',  sim_aParm4)

  ! chemistry
  call RuntimeParameters_get('sim_tdust', sim_tdust)
  call RuntimeParameters_get('sim_init_Hp', sim_init_Hp)

  ! stratbox, read in kolmogorov velocity data
  call RuntimeParameters_get('sim_velcubeFile', sim_velcubeFile)
  call RuntimeParameters_get('sim_machTurb', sim_machTurb)
  call RuntimeParameters_get('xmax', sim_xMax) ! declared by Grid/GridMain
  call RuntimeParameters_get('xmin', sim_xMin) ! needed to do grid interp
  call RuntimeParameters_get('ymax', sim_yMax)
  call RuntimeParameters_get('ymin', sim_yMin)
  call RuntimeParameters_get('zmax', sim_zMax)
  call RuntimeParameters_get('zmin', sim_zMin)

  ! for sim_cubeFile read and broadcast
  call Driver_getMype(GLOBAL_COMM, sim_MyPE)
  call Driver_getComm(GLOBAL_COMM, sim_Comm)

  ! read number size of the cube file and communicate it
  if (sim_myPE .eq. MASTER_PE) then
    open(unit = 55, file = sim_velcubeFile, status = 'old')
    read(55,*) dumstr, sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)
  endif
  call MPI_Bcast(sim_nCD, MDIM, MPI_INTEGER, MASTER_PE, sim_comm, istat)

  ! allocate field arrays for cubeFile data
  allocate(sim_velxArr(sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_velxArr")
  allocate(sim_velyArr(sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_velyArr")
  allocate(sim_velzArr(sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_velzArr")

  if (sim_myPE .eq. MASTER_PE) then
    do i = 1,sim_nCD(IAXIS)
      do j = 1,sim_nCD(JAXIS)
        do k = 1,sim_nCD(KAXIS)
          read(55,*) ii,jj,kk, sim_velxArr(i,j,k), sim_velyArr(i,j,k), &
                     sim_velzArr(i,j,k)
        enddo
      enddo
    enddo
    close(55)

    ! ----------------------
    ! set velocity magnitude
    ! just an experiment,
    ! do not use in production runs!!!!!
    ! ----------------------
    ! the gas layer has uniform temperature and hence uniform internal energy,
    ! because P and rho follow same z-profile
    eint = sim_p/sim_rho/(5./3-1.)  ! sim_gamma not yet initialized because we are using Multispecies

    ! mean specific KE
    ek = sum(sim_velxArr(:,:,:)**2+sim_velyArr(:,:,:)**2+sim_velzArr (:,:,:)**2) &
         / (sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS))

    print *, "[Simulation_init] WARNING: using approx velocity init with sim_Machturb"
    print *, "[Simulation_init]    various prefactors of order unity not accounted for,"
    print *, "[Simulation_init]    sim_gamma=5/3 hardcoded in."
    print *, "[Simulation_init]    incorrect turb velocity for external IGM (but that shouldn't matter)"
    print *, "[Simulation_init]    DO NOT USE THESE RUNS FOR PUBLICATION-QUALITY SCIENCE!!"
    print *, "[Simulation_init]    - AT, 20190321"

    print *, "[Simulation_init] sim_machTurb", sim_machTurb
    print *, "[Simulation_init] eint",eint, "ek",ek, "before renorm"

    ! turbulent mach number is ~ <v> / <c_s> ~ sqrt(ek/eint)
    ! I am being VERY lazy about factors of unity
    ! mainly, eint is not <c_s^2> in general
    ! but OK for our purposes, we care about energy ratio.
    ! ALSO: for the stratbox setup, this sets wrong velocity for IGM (will be
    ! too low compared to what is expected...although IDK how much we know about
    ! IGM velocity anyways)
    sim_velxArr = sim_velxArr * sim_machTurb / (ek/eint)**0.5  ! Mach_target / Mach_current / sqrt(3)
    sim_velyArr = sim_velyArr * sim_machTurb / (ek/eint)**0.5
    sim_velzArr = sim_velzArr * sim_machTurb / (ek/eint)**0.5

    ! recompute for debugging purposes
    ek = sum(sim_velxArr(:,:,:)**2+sim_velyArr(:,:,:)**2+sim_velzArr (:,:,:)**2) &
         / (sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS))
    print *, "[Simulation_init] eint",eint, "ek",ek, "after renorm"
    print *, "[Simulation_init] eint**0.5 in km/s =", eint**0.5 / 1e5
    print *, "[Simulation_init] ek**0.5 in km/s =", ek**0.5 / 1e5
    print *, "[Simulation_init]"
    print *, "[Simulation_init] ek/eint", ek/eint
    print *, "[Simulation_init] ek**0.5 / eint**0.5 =", (ek/eint)**0.5
    !call Driver_abortFlash("Done testing!")

  endif

  call MPI_Bcast(sim_velxArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_velyArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_velzArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)

end subroutine Simulation_init

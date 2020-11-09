!!****if* source/Simulation/SimulationMain/Cube/Simulation_init
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
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Also initializes initial conditions for cube problem
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getComm

  implicit none
#include "Flash.h"
#include "Eos.h"
#include "constants.h"
#include "Flash_mpi.h"


  integer :: i,j,k, istat, ii,jj,kk
  character(len=255),save :: dumstr

  call Driver_getMype(GLOBAL_COMM, sim_MyPE)
  call Driver_getComm(GLOBAL_COMM, sim_Comm)


  call RuntimeParameters_get( 'sim_cubeFile' , sim_cubeFile )

  call RuntimeParameters_get( 'gamma', sim_gamma)
  call RuntimeParameters_get( 'smlrho', smlrho)
  call RuntimeParameters_get( 'smallp', smallp)
  call RuntimeParameters_get( 'smallX', smallX)
  call RuntimeParameters_get( 'xmax', sim_xMax)
  call RuntimeParameters_get( 'xmin', sim_xMin)
  call RuntimeParameters_get( 'ymax', sim_yMax)
  call RuntimeParameters_get( 'ymin', sim_yMin)
  call RuntimeParameters_get( 'zmax', sim_zMax)
  call RuntimeParameters_get( 'zmin', sim_zMin)
  

  call PhysicalConstants_get( 'Boltzmann', sim_boltz)
  call PhysicalConstants_get( 'proton mass', sim_mH)
  call PhysicalConstants_get( 'pi', sim_pi)

  call RuntimeParameters_get('eos_singleSpeciesA', sim_molarMass)
  call PhysicalConstants_get('ideal gas constant', sim_gasconstant)
  call PhysicalConstants_get('proton mass', sim_protonMass)

  call RuntimeParameters_get('sim_init_Hp', sim_init_Hp)
  
  call RuntimeParameters_get( 'bx0', sim_magx)
  call RuntimeParameters_get( 'by0', sim_magy)
  call RuntimeParameters_get( 'bz0', sim_magz)

  call RuntimeParameters_get( 'sim_tdust', sim_tdust)
  call RuntimeParameters_get('killdivb',sim_killdivb)
! New Parameterized heating and cooling
!  call RuntimeParameters_get('sim_Z',  sim_Z)
!  call RuntimeParameters_get('sim_G0', sim_G0)
!  call RuntimeParameters_get('h_uv', sim_pe_h)
!  call RuntimeParameters_get('h_cr', sim_cr_h)
!  call RuntimeParameters_get('use_constant_heating', sim_constant_heating)
!  call RuntimeParameters_get('stratifyHeat',   sim_stratify_heating)

  sim_abar = 1.0 + sim_abundM*sim_metal


  ! read number size of the cube file and communicate it
  if (sim_myPE .eq. MASTER_PE) then
    open(unit = 55, file = sim_cubeFile, status = 'old')
    read(55,*) dumstr, sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)
  endif
  call MPI_Bcast(sim_nCD, MDIM, MPI_INTEGER, MASTER_PE, sim_comm, istat)

  ! allocate field arrays for cubeFile data
  allocate(sim_densArr(sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_densArr")
  allocate(sim_presArr(sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_presArr")
  allocate(sim_gpotArr(sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_gpotArr")
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
          read(55,*) ii,jj,kk, sim_densArr(i,j,k), sim_presArr(i,j,k), &
          & sim_velxArr(i,j,k), sim_velyArr(i,j,k), sim_velzArr(i,j,k), &
          & sim_gpotArr(i,j,k)
          !print *, i,ii, j, jj, k, kk, sim_densArr(i,j,k), sim_presArr(i,j,k), sim_gpotArr(i,j,k)
        enddo
!        read(55,*)
      enddo
!      read(55,*)
    enddo
    close(55)
  endif

  call MPI_Bcast(sim_densArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_presArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_gpotArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_velxArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_velyArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_velzArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS), &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)


end subroutine Simulation_init

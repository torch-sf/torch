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
  use Logfile_interface, ONLY : Logfile_stamp
  
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
  call RuntimeParameters_get('killdivb', sim_killdivb)
! New Parameterized heating and cooling
!  call RuntimeParameters_get('sim_Z',  sim_Z)
!  call RuntimeParameters_get('sim_G0', sim_G0)
!  call RuntimeParameters_get('h_uv', sim_pe_h)
!  call RuntimeParameters_get('h_cr', sim_cr_h)
!  call RuntimeParameters_get('use_constant_heating', sim_constant_heating)
!  call RuntimeParameters_get('stratifyHeat',   sim_stratify_heating)

  ! stratbox, gravity profile parameters, stellar disk
  call RuntimeParameters_get('sim_aParm1',  sim_aParm1)
  call RuntimeParameters_get('sim_aParm2',  sim_aParm2)
  call RuntimeParameters_get('sim_aParm3',  sim_aParm3)
  call RuntimeParameters_get('sim_aParm4',  sim_aParm4)
  call RuntimeParameters_get('sim_withStaticGrav', sim_withStaticGrav)

  ! VorAMR switches and info
  call RuntimeParameters_get('use_voramr', use_voramr)
  call RuntimeParameters_get('voramr_source', voramr_source)
  call RuntimeParameters_get('voramr_input', voramr_input)
  call RuntimeParameters_get('use_localRef', use_localRef)
  call RuntimeParameters_get('center_localRef', center_localRef)
  call RuntimeParameters_get('localRef_x', localRef_x)
  call RuntimeParameters_get('localRef_y', localRef_y)
  call RuntimeParameters_get('localRef_z', localRef_z)
  call RuntimeParameters_get('localRef_r', localRef_r)
  call RuntimeParameters_get('refine_on_particle_count', refPartCount)
  if (sim_myPE.EQ.MASTER_PE) then
     if (refPartCount) call Logfile_stamp('refine_on_particle_count is set to .TRUE.', "[Simulation_init]")
     if (sim_withStaticGrav) call Logfile_stamp('sim_withStaticGrav is set to .TRUE.', "[Simulation_init]")
     if (use_voramr) then
        call Logfile_stamp('use_voramr is set to .TRUE. proceeding with VorAMR.', "[Simulation_init]")
        call Logfile_stamp(voramr_source, "[Simulation_init]")
     endif
     if (use_localRef) call Logfile_stamp('use_localRef is set to .TRUE. only refining part of grid initially', &
          & "[Simulation_init]")
  endif
  
   !! Derefinement outside rectangular region of interest
  call RuntimeParameters_get('use_deref', use_deref)
  call RuntimeParameters_get('deref_lref', deref_lref)
  call RuntimeParameters_get('deref_xl', deref_xl)
  call RuntimeParameters_get('deref_xr', deref_xr)
  call RuntimeParameters_get('deref_yl', deref_yl)
  call RuntimeParameters_get('deref_yr', deref_yr)
  call RuntimeParameters_get('deref_zl', deref_zl)
  call RuntimeParameters_get('deref_zr', deref_zr)
  if (sim_myPE.EQ.MASTER_PE) then
     if (use_deref) call Logfile_stamp('derefining outside region of interest', "[Simulation_init]")
  endif
  
#ifdef ELEMENTS
  call RuntimeParameters_get('nelements', sim_nelements)
  
! Check for consistent setup.
  if (NMASS_SCALARS .ne. sim_nelements) then
    if (sim_myPE .eq. MASTER_PE) then
      print*, "Number of mass scalars is not consitent with sim_nelements: NMASS_SCALARS, sim_nelements = ", NMASS_SCALARS, sim_nelements
    end if
    call Driver_abortFlash("Inconsitency between number of elements and mass scalars. Check flash.par and Simulation Config.")
  endif
#endif

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
  
#ifdef ELEMENTS
  allocate(sim_elemArr(NMASS_SCALARS, sim_nCD(IAXIS), sim_nCD(JAXIS), sim_nCD(KAXIS)), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_elemArr")
#endif

  if (sim_myPE .eq. MASTER_PE) then
    do i = 1,sim_nCD(IAXIS)
      do j = 1,sim_nCD(JAXIS)
        do k = 1,sim_nCD(KAXIS)
#ifdef ELEMENTS
          read(55,*) ii,jj,kk, sim_densArr(i,j,k), sim_presArr(i,j,k), &
          & sim_velxArr(i,j,k), sim_velyArr(i,j,k), sim_velzArr(i,j,k), &
          & sim_gpotArr(i,j,k), sim_elemArr(:,i,j,k)
#else
          read(55,*) ii,jj,kk, sim_densArr(i,j,k), sim_presArr(i,j,k), &
          & sim_velxArr(i,j,k), sim_velyArr(i,j,k), sim_velzArr(i,j,k), &
          & sim_gpotArr(i,j,k)
          !print *, i,ii, j, jj, k, kk, sim_densArr(i,j,k), sim_presArr(i,j,k), sim_gpotArr(i,j,k)
#endif
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
#ifdef ELEMENTS
  call MPI_Bcast(sim_elemArr, sim_nCD(IAXIS)*sim_nCD(JAXIS)*sim_nCD(KAXIS)*NMASS_SCALARS, &
  & FLASH_REAL, MASTER_PE, sim_comm, istat)
#endif


end subroutine Simulation_init

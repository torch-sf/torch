!!****if* source/Particles/ParticlesInitialization/Particles_init
!!
!! NAME
!!    Particles_init
!!
!! SYNOPSIS
!!    Particles_init( logical(in) :: restart )
!!
!! DESCRIPTION
!!
!!    General initialization routine for the particle module.
!!
!! ARGUMENTS
!!
!!    restart  : indicates if run is starting from scratch or restarting
!!               from checkpoint file
!!
!! PARAMETERS
!!
!!    useParticles   BOOLEAN [TRUE]  Should particles be used in this simulation?
!!    pt_maxPerProc  INTEGER [100]   Maximum number of particles per processor. Allocates array space
!!                                   Particles are distributed per PROCESSOR rather than per BLOCK
!!    pt_dtFactor    REAL    [0.5]   Factor to make sure that time step is small enough that particles
!!    pt_dtChangeTolerance REAL [0.4] For uncorrected Estimated Midpoint propagation scheme:
!!                                    Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities
!!    pt_small       REAL    [1.0E-10] Used for general comparisons of real values 
!!                                   For example, IF (abs(real1 - real2) .lt. pt_small) THEN
!!                                   don't move farther than one block in each step
!!    
!!  NOTE 
!!
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_init ( restart)
  
  use Particles_data

  use Driver_interface, ONLY : Driver_abortFlash, Driver_getComm, Driver_getMype,&
       Driver_getNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
    RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use pt_interface, ONLY : pt_mapStringParamToInt, pt_picInit, pt_setMask, pt_setDataStructures
  use Particles_interface, ONLY : Particles_specifyMethods, Particles_sinkInit, Particles_rayInit
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "GridParticles.h"
  logical, INTENT(in) :: restart

  integer :: ierr
  integer :: i,j
  logical :: isLattice, isWithDensity

!-------------------------------------------------------------------------------
  

  call RuntimeParameters_get ( "useParticles", useParticles)
  call RuntimeParameters_get ( "keepLostParticles", pt_keepLostParticles)

  !! get out of here if we don't need Particles

  !! but the following initializations are neccessary if the Particles
  !! unit is included in the simulation for Particles_sendOutputData 
  !(called from IO_setScalar) to get the log message correct.

  pt_numLocal = 0   ! number of particles on this processor
  call Driver_getMype(MESH_COMM, pt_meshMe)
  call Driver_getNumProcs(MESH_COMM, pt_meshNumProcs)
  call Driver_getComm(MESH_COMM, pt_meshComm)

  call Driver_getMype(GLOBAL_COMM, pt_globalMe)
  call Driver_getNumProcs(GLOBAL_COMM, pt_globalNumProcs)
  call Driver_getComm(GLOBAL_COMM, pt_globalComm)

  pt_restart = restart   ! Is this a restart?

     
  if (.NOT. useParticles) return

  
  ! Time stepping parameters, for Particles_computeDt
  call RuntimeParameters_get("pt_dtFactor", pt_dtFactor)
  
  ! Location bounds, needed for Particles_initPosition and Particles_mapParticles2Mesh
  call RuntimeParameters_get( "xmin", pt_xmin)
  call RuntimeParameters_get( "xmax", pt_xmax)
  call RuntimeParameters_get( "ymin", pt_ymin)
  call RuntimeParameters_get( "ymax", pt_ymax)
  call RuntimeParameters_get( "zmin", pt_zmin)
  call RuntimeParameters_get( "zmax", pt_zmax)

    ! space allocation, for Particles_init
  call RuntimeParameters_get ("pt_maxPerProc", pt_maxPerProc)

  ! geometry, for Particles_initPositions
  call RuntimeParameters_get("geometry",pt_str_geometry)
  call RuntimeParameters_mapStrToInt(pt_str_geometry, pt_geometry)

  ! defined variously for the advance methods, used in Particles_initAttributes

  !  Set a small comparison number, for timestep comparisons in some Particles_advance
  call RuntimeParameters_get ("pt_small", pt_small)
  call RuntimeParameters_get ("pt_dtChangeTolerance", pt_dtChangeTolerance)
  call RuntimeParameters_get ("pt_numAtOnce", pt_numAtOnce)
  call RuntimeParameters_get ("pt_logLevel", pt_logLevel)
  call RuntimeParameters_get ("pt_resetTag", pt_resetTag)


  ! We create Particles_specifyMethods() at setup time, which is used 
  ! to populate PART_TYPE, PART_INITMETHOD, PART_MAPMETHOD in pt_typeInfo.
  pt_typeInfo=0
  call Particles_specifyMethods()
  pt_typeInfo(PART_LOCAL_MAX,:)=pt_maxPerProc
  pt_numLost=0

  !! If any of the particle types are initialized with the lattice method
  !! then the following parameters are needed.
  isLattice=.false.
  isWithDensity=.false.

  do i=1,NPART_TYPES
     isLattice=isLattice.or.(pt_typeInfo(PART_INITMETHOD,i)==LATTICE)
     isWithDensity=isWithDensity.or.&
          (pt_typeInfo(PART_INITMETHOD,i)==CELLMASS).or.&
          (pt_typeInfo(PART_INITMETHOD,i)==REJECTION).or.&
          (pt_typeInfo(PART_INITMETHOD,i)==WITH_DENSITY)
  end do

  if(isLattice) then
#ifndef PART_INITMETHOD_LATTICE
     call Driver_abortFlash("Particles are set up to use a Lattice initialization method,"// &
          " but a unit to implement that method has not been compiled in!")
#endif
     call RuntimeParameters_get ("pt_numX", pt_numX)
     call RuntimeParameters_get ("pt_numY", pt_numY)
     call RuntimeParameters_get ("pt_numZ", pt_numZ)
     call RuntimeParameters_get ("pt_initialXMin", pt_initialXMin)
     call RuntimeParameters_get ("pt_initialXMax", pt_initialXMax)
     call RuntimeParameters_get ("pt_initialYMin", pt_initialYMin)
     call RuntimeParameters_get ("pt_initialYMax", pt_initialYMax)
     call RuntimeParameters_get ("pt_initialZMin", pt_initialZMin)
     call RuntimeParameters_get ("pt_initialZMax", pt_initialZMax)
     call RuntimeParameters_get ("pt_initialRadius", pt_initialRadius)
  end if

  if(isWithDensity) then
#ifndef PART_INITMETHOD_WITHDENSITY
     call Driver_abortFlash("Particles are set up to use a density-based initialization method,"// &
          " but a unit to implement that method has not been compiled in!")
#endif
     call RuntimeParameters_get ("pt_pRand", pt_pRand)
     call RuntimeParameters_get("pt_numParticlesWanted", pt_numParticlesWanted)
  end if

    ! Reduce guard cell fills
  call RuntimeParameters_get ("reduceGcellFills", pt_reduceGcellFills)
#if NPART_TYPES != 1 || !defined (PASSIVE_PART_TYPE)
  if (pt_reduceGcellFills) then !For now, works with passive particles only
     call Logfile_stamp(&
&'Particle types other than PASSIVE are included in the configuration, therefore ignoring reduceGcellFills for Particles_advance!'&
          ,'[Particles_init]')
     pt_reduceGcellFills = .false.
  end if
#endif

  
  !! Call the routine that can set the mask for guardcell fill when asked by the Particles unit
  call pt_setMask()
  

  if (.not. restart) then !if we starting from scratch
     
     ! allocate data structures for particle storage
     if (.not. allocated(particles)) &
          allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
     if (ierr /= 0) then
        call Driver_abortFlash("Particles_init:  could not allocate particle array")
     endif
          
     !initialize all particles to NONEXISTENT
     particles = NONEXISTENT
     
  end if  ! end of .not. restart

  
  ! Added so I can do a few operations (like initialize the local tag
  ! numbers for each proc properly when keeping the same tags for
  ! both sinks and particles) in simulations. - JW
  ! Note this is always initialized to NONEXISTANT.
  ! Call pt_gatherGlobal to properly fill.
  if (.not. allocated(allproc_particles)) &
       allocate (allproc_particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) then
     call Driver_abortFlash("Particles_init:  could not allocate particle array")
  endif
  allproc_particles = NONEXISTENT

  ! allocate data structures for AMUSE to pick up new particles. - JW
  if (.not. allocated(new_massive_tags)) &
      allocate (new_massive_tags(pt_maxPerProc), stat=ierr)
  if (ierr /= 0) then
    call Driver_abortFlash("Particles_init:  could not allocate new_tags")
  endif
 
     number_new_massive = 0

  call pt_setDataStructures()
  

  !! The index list is created for interaction with the GridParticles subunit.
  !! The GridParticles subunit is designed to move entities between blocks when
  !! those entities are modeling some physics that does not fit with the Eulerian
  !! data structures native to the grid, but need to interact with one or more
  !! grid variables. The subunit also maps attributes associated with these
  !! entities to and from the mesh. In order to allow GridParticles to be agnostic
  !! to the attributes of the various types of entities (particles vs rays for example)
  !! the client code unit initializes the index map; this section does that for
  !! Lagrangian particles modeled by the particles unit.

  pt_indexList(1:GRPT_ALL)=GRPT_RESET
  pt_indexCount=GRPT_ALL

  pt_indexList(GRPT_POSX_IND)=POSX_PART_PROP
  pt_indexList(GRPT_POSY_IND)=POSY_PART_PROP
  pt_indexList(GRPT_POSZ_IND)=POSZ_PART_PROP
  pt_indexList(GRPT_BLK_IND) =BLK_PART_PROP
  pt_indexList(GRPT_TAG_IND) =TAG_PART_PROP
  pt_indexList(GRPT_PROC_IND)=PROC_PART_PROP
#ifdef PART_VIRTUAL
  pt_indexList(GRPT_VIRTUAL)=VIRT_PART_PROP
#endif

#ifdef TMPX_PART_PROP  
  !! If any other position co-ordinates exist for holding
  !! temporary values that need to have boundary conditions
  !! applied to them
#ifdef POSPREDX_PART_PROP
  call Driver_abortFlash(&
       "Particles_init: two types of temporary position attributes are not supported yet (TMPX_PART_PROP / POSPREDX_PART_PROP)")
#endif
  pt_indexList(GRPT_POSXTMP_IND)=TMPX_PART_PROP
  pt_indexList(GRPT_POSYTMP_IND)=TMPY_PART_PROP
  pt_indexList(GRPT_POSZTMP_IND)=TMPZ_PART_PROP
  pt_indexList(GRPT_POSTMP) = GRPT_EXIST
#endif

#ifdef POSPREDX_PART_PROP
  pt_indexList(GRPT_POSXTMP_IND)=POSPREDX_PART_PROP
  pt_indexList(GRPT_POSYTMP_IND)=POSPREDY_PART_PROP
  pt_indexList(GRPT_POSZTMP_IND)=POSPREDZ_PART_PROP
  pt_indexList(GRPT_POSTMP) = GRPT_EXIST
#endif

#ifdef VELX_PART_PROP
  pt_indexList(GRPT_VELX_IND)=VELX_PART_PROP
  pt_indexList(GRPT_VELY_IND)=VELY_PART_PROP
  pt_indexList(GRPT_VELZ_IND)=VELZ_PART_PROP
  pt_indexList(GRPT_VEL) = GRPT_EXIST
#endif

#ifdef VELPREDX_PART_PROP
  pt_indexList(GRPT_VELXTMP_IND)=VELX_PART_PROP
  pt_indexList(GRPT_VELYTMP_IND)=VELY_PART_PROP
  pt_indexList(GRPT_VELZTMP_IND)=VELZ_PART_PROP
  pt_indexList(GRPT_VELTMP) = GRPT_EXIST
#endif


  call pt_picInit()
  ! initialize sink particles
  call Particles_sinkInit(restart)
  ! initialize ray tracer
  call Particles_rayInit(restart)
  return

end subroutine Particles_init

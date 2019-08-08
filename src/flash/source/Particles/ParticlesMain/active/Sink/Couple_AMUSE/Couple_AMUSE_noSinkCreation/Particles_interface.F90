!!****h* source/Particles/Particles_interface
!!
!! This is the header file for the particles module that defines its
!! public interfaces.
!!***
Module Particles_interface

  implicit none

#include "Flash.h"
#include "constants.h"

  interface
    subroutine Particles_advance (dtOld,dtNew)
      real, INTENT(in) :: dtOld, dtNew
    end subroutine Particles_advance
  end interface

  interface
    subroutine Particles_computeDt (blockID, dt_part, dt_minloc)
      real, INTENT(inout)    :: dt_part
      integer, INTENT(inout) :: dt_minloc(5)
      integer, INTENT(in)    :: blockID
    end subroutine Particles_computeDt
  end interface

  interface
    subroutine Particles_dump(blockCount,blockList,nstep,time,dt)
      integer, intent(IN) :: blockCount
      integer, intent(IN) :: blockList(blockCount)
      integer, intent(IN) :: nstep
      real, intent(IN)    :: time, dt
    end subroutine Particles_dump
  end interface

  interface
    subroutine Particles_finalize()
    end subroutine Particles_finalize
  end interface

  interface
    subroutine Particles_getGlobalNum(globalNumParticles)
      integer, intent(out)  :: globalNumParticles
    end subroutine Particles_getGlobalNum
  end interface

  interface
    subroutine Particles_getLocalNum(blockID, localNumParticles)
      integer, intent(in) :: blockID
      integer, intent(out)  :: localNumParticles
    end subroutine Particles_getLocalNum
  end interface

  interface
    subroutine Particles_init( restart)
      logical, intent(IN) :: restart
    end subroutine Particles_init
  end interface

  interface
    subroutine Particles_initPositions ( partPosInitialized,updateRefine)
      logical, INTENT(INOUT) :: partPosInitialized
      logical, INTENT(OUT) :: updateRefine
    end subroutine Particles_initPositions
  end interface

  interface
    subroutine Particles_initData(restart,partPosInitialized)
      logical, INTENT(IN) :: restart
      logical,INTENT(INOUT) :: partPosInitialized
    end subroutine Particles_initData
  end interface

  interface
    subroutine Particles_longRangeForce (particles,p_count,mapType,pot_var)

      integer, intent(IN) :: p_count,mapType
      real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles
      integer, intent(IN), optional :: pot_var
    end subroutine Particles_longRangeForce
  end interface

  interface
    subroutine Particles_putLocalNum(localNumParticles)
      integer, intent(in)  :: localNumParticles
    end subroutine Particles_putLocalNum
  end interface

  interface
    subroutine Particles_sendOutputData()
    end subroutine Particles_sendOutputData
  end interface

  interface
    subroutine Particles_shortRangeForce ()
    end subroutine Particles_shortRangeForce
  end interface

  interface
     subroutine Particles_unitTest(fileUnit,perfect)
       integer, intent(in) ::  fileUnit
       logical, intent(INOUT) :: perfect
     end subroutine Particles_unitTest
  end interface

  interface
     subroutine Particles_updateRefinement(oldLocalNumBlocks)
       integer,intent(INOUT) :: oldLocalNumBlocks
     end subroutine Particles_updateRefinement
  end interface

  interface
     subroutine Particles_mapFromMesh (mapType,numAttrib, attrib, pos, bndBox,&
          deltaCell,solnVec,  partAttribVec)
       
       integer, INTENT(in) :: mapType,numAttrib
       integer, dimension(2, numAttrib),intent(IN) :: attrib
       real,dimension(MDIM), INTENT(in)    :: pos,deltaCell
       real, dimension(LOW:HIGH,MDIM), intent(IN) :: bndBox
       real, pointer       :: solnVec(:,:,:,:)
       real,dimension(numAttrib), intent(OUT) :: partAttribVec
       
     end subroutine Particles_mapFromMesh
  end interface
  
  interface
     subroutine Particles_mapToMeshOneBlk(blkLimitsGC,guard, blockID,&
          particles,numParticles,pt_attribute,buff)
       implicit none
       integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: blkLimitsGC
       integer, dimension(MDIM),intent(IN) :: guard
       integer, intent(IN) :: blockID
       integer, intent(in) :: numParticles
       real, intent(in) :: particles(NPART_PROPS,numParticles)
       integer, intent(IN) :: pt_attribute
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
            INTENT(INOUT) :: buff
     end subroutine Particles_mapToMeshOneBlk
  end interface

  interface
     subroutine Particles_updateAttributes()
       implicit none
     end subroutine Particles_updateAttributes
  end interface

  interface
     subroutine Particles_updateGridVar(partProp, varGrid, mode)
       integer, INTENT(in) :: partProp, varGrid
       integer, INTENT(in), optional :: mode
     end subroutine Particles_updateGridVar
  end interface

  interface
     subroutine Particles_accumCount(var)
       integer,intent(IN) :: var
     end subroutine Particles_accumCount
  end interface

  interface
     subroutine Particles_getCountPerBlk(perBlkCount)
       
       integer,dimension(MAXBLOCKS),intent(OUT) :: perBlkCount
     end subroutine Particles_getCountPerBlk
  end interface

  interface
     subroutine Particles_initForces()
     end subroutine Particles_initForces
  end interface

  interface
     subroutine Particles_specifyMethods()
     end subroutine Particles_specifyMethods
  end interface

  interface
     subroutine Particles_manageLost(mode)
       integer, intent(IN) :: mode
     end subroutine Particles_manageLost
  end interface

  interface
      subroutine Particles_clean()
      end subroutine Particles_clean
   end interface

  interface
     subroutine Particles_addNew (count, pos, success)
       integer, INTENT(in) :: count
       real, optional, dimension(MDIM,count), intent(IN)::pos
       logical, intent(OUT) :: success
     end subroutine Particles_addNew
  end interface
  
  interface
    subroutine Particles_sinkAccelGasOnSinks(accelProps)
      implicit none
      integer, intent(in), OPTIONAL :: accelProps(MDIM)
    end subroutine Particles_sinkAccelGasOnSinks
  end interface
  
  interface
    subroutine Particles_sinkAccelSinksOnGas(blockCount,blockList, accelVars)
      implicit none
      integer,intent(IN) :: blockCount
      integer,dimension(blockCount),intent(IN) :: blockList
      integer, intent(in), OPTIONAL :: accelVars(MDIM)
    end subroutine Particles_sinkAccelSinksOnGas
  end interface

  interface
    subroutine Particles_sinkAdvanceParticles(dr_dt)
      real, INTENT(in) :: dr_dt
    end subroutine Particles_sinkAdvanceParticles
  end interface
  
  interface
    subroutine Particles_sinkComputeDt(blockID,dt_sink,dt_minloc)
      integer, INTENT(in)    :: blockID
      real, INTENT(inout)    :: dt_sink
      integer, INTENT(inout) :: dt_minloc(5)
    end subroutine Particles_sinkComputeDt
  end interface
  
  interface
    subroutine Particles_sinkCreateAccrete(dt)
      real, intent(IN) :: dt
    end subroutine Particles_sinkCreateAccrete
  end interface 
  
  interface
    subroutine Particles_sinkInit(restart)
      logical, intent(in) :: restart
    end subroutine Particles_sinkInit
  end interface
  
  interface
    subroutine Particles_sinkMarkRefineDerefine()
    end subroutine Particles_sinkMarkRefineDerefine
  end interface
  
  interface
    subroutine Particles_sinkMoveParticles(regrid)
      logical, intent(in) :: regrid
    end subroutine Particles_sinkMoveParticles
  end interface
  
  interface
    subroutine Particles_sinkSortParticles()
    end subroutine Particles_sinkSortParticles
  end interface
  
  interface
    subroutine Particles_sinkSumAttributes(sums, attribs, factor)
      implicit none
      real,intent(OUT)   :: sums(:)
      integer,intent(in) :: attribs(:)
      integer,intent(in),OPTIONAL :: factor
    end subroutine Particles_sinkSumAttributes
  end interface
  
  interface
    subroutine Particles_sinkSyncWithParticles(sink_to_part)
      logical, intent(in) :: sink_to_part
    end subroutine Particles_sinkSyncWithParticles
  end interface

! ray tracing interfaces
  interface
    subroutine Particles_rayAdvance(dt)
      implicit none
      real, intent(IN) :: dt
    end subroutine Particles_rayAdvance
  end interface

  interface
    subroutine Particles_rayInit(restart)
      implicit none
      logical, intent(in) :: restart
    end subroutine Particles_rayInit
  end interface

! stars particles interfaces
  interface
    subroutine Particles_starsAdvance(dt)
      implicit none
      real, intent(IN) :: dt
    end subroutine Particles_starsAdvance
  end interface

  interface
    subroutine Particles_starsInit(restart)
      implicit none
      logical, intent(in) :: restart
    end subroutine Particles_starsInit
  end interface
   
!!! AMUSE additions. - JW

  interface
    subroutine Particles_MassRemoval(dt, mass_acc, xloc, yloc, zloc, velx, vely, velz)
      real, intent(IN) :: dt
	  real, dimension(10), intent(INOUT) :: mass_acc, velx, vely, velz, xloc, yloc, zloc
    end subroutine Particles_MassRemoval
  end interface

  interface
    subroutine dens_removal(made_stars)
      implicit none
      logical, intent(out) :: made_stars
    end subroutine dens_removal
  end interface

!  interface
!    subroutine gaussian(x, gaus, mu_in, var_in)
!      #include "constants.h"
!      implicit none
!     real, dimension(:), intent(in) :: x
!     real, optional, intent(in) :: mu_in, var_in ! Center of gaus, variance
!     real, dimension(:), allocatable, intent(out) :: gaus
!    end subroutine gaussian
!  end interface

  interface
    subroutine Particles_sinkKickGas(blockCount,blockList, dt, force_sum_x_glob, force_sum_y_glob, force_sum_z_glob, accelVars)
      implicit none
      real, intent(IN) :: dt
      integer,intent(IN) :: blockCount
      integer,dimension(blockCount),intent(IN) :: blockList
      real, intent(INOUT) :: force_sum_x_glob, force_sum_y_glob, force_sum_z_glob
      integer, intent(in), OPTIONAL :: accelVars(MDIM)
    end subroutine Particles_sinkKickGas
  end interface

  interface
    subroutine Particles_moveAndSort(regrid)
      implicit none
      logical, intent(in) :: regrid
    end subroutine Particles_moveAndSort
  end interface
  
! A method to inject kinetic and thermal energy into the grid. - JW
  interface
    subroutine Particles_energyInjection(energy, fracKin, mass, xloc, yloc, zloc, dt)
 
      implicit none
      real, intent(in) :: energy, mass, xloc, yloc, zloc
      real, intent(inout) :: fracKin, dt
    end subroutine Particles_energyInjection
  end interface

  interface
    subroutine Particles_sinkMakeStars(dt)
      real, intent(IN) :: dt
    end subroutine Particles_sinkMakeStars
  end interface

  interface
    subroutine Particles_sinkWind(dt)
    
      implicit none
      real, intent(in) :: dt

    end subroutine Particles_sinkWind
  end interface

  interface
    subroutine Particles_wind(starMass, mass, xloc, yloc, zloc, dt)
 
      implicit none
      real, intent(in) :: starMass, mass, xloc, yloc, zloc
      real, intent(inout) ::  dt
    end subroutine Particles_wind
  end interface

end Module Particles_interface

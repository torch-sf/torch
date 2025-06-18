!!****h* source/Simulation/Simulation_interface
!!
!! This is the header file for the Simulation module
!! that defines its public interfaces.
!!***
Module Simulation_interface
  implicit none
#include "constants.h"
  interface
     subroutine Simulation_defineDomain(initialDomain,boundaries,nblks)
       implicit none
       integer,dimension(MDIM),intent(IN) :: nblks
       integer,dimension(2*MDIM,nblks(IAXIS),nblks(JAXIS),nblks(KAXIS)),&
            intent(OUT)::boundaries
       logical,dimension(nblks(IAXIS),nblks(JAXIS),nblks(KAXIS)),&
            intent(OUT)::initialDomain
     end subroutine Simulation_defineDomain
  end interface

  interface
     subroutine Simulation_finalize()
       implicit none
     end subroutine Simulation_finalize
  end interface

  interface
     subroutine Simulation_getRenormGroup(mscalar,group)
       implicit none
       integer, intent(out) ::group
       integer, intent(in) :: mscalar
     end subroutine Simulation_getRenormGroup
  end interface

  interface
     PURE subroutine Simulation_getVarnameType(varname,vartype)
       implicit none
       integer, intent(out) :: vartype
       integer, intent(in) :: varname
     end subroutine Simulation_getVarnameType
  end interface

  interface
     subroutine Simulation_initBlock(blockID)
       implicit none
       integer, intent(in) :: blockID
     end subroutine Simulation_initBlock
  end interface

  interface
     subroutine Simulation_init()
       implicit none
     end subroutine Simulation_init
  end interface

  interface
     subroutine Simulation_initParticleAttrib(restart)
       logical,intent(in) :: restart
     end subroutine Simulation_initParticleAttrib
  end interface

  interface
     subroutine Simulation_initSpecies()
       implicit none
     end subroutine Simulation_initSpecies
  end interface

  interface
     subroutine Simulation_mapIntToStr(key, str, block)
       implicit none
       integer, intent(in) :: key, block
       character(len=*), intent(inout) :: str
     end subroutine Simulation_mapIntToStr
  end interface

  interface
     subroutine Simulation_mapStrToInt(str,key,map)
       implicit none
       character(len=*), intent(in) :: str
       integer, intent(out) :: key
       integer, intent(in) :: map
     end subroutine Simulation_mapStrToInt
  end interface

  interface
     subroutine Simulation_sendOutputData()
       implicit none
     end subroutine Simulation_sendOutputData
  end interface

  interface
     subroutine Simulation_mapParticlesVar(part_key, var_key, var_type)
       implicit none
       integer, intent(in)  :: part_key
       integer, intent(out) :: var_key, var_type

     end subroutine Simulation_mapParticlesVar
  end interface

  interface
     subroutine Simulation_initRestart()
       implicit none

     end subroutine Simulation_initRestart
  end interface

  interface
     subroutine Simulation_customizeProlong(beforeOrAfter)
       implicit none
       integer, intent (IN) :: beforeOrAfter

     end subroutine Simulation_customizeProlong
  end interface

  interface
     subroutine Simulation_computeAnalytical(blockID,  tcurr)
      implicit none
      integer, intent (IN) :: blockID
      real   , intent (IN) :: tcurr
     end subroutine Simulation_computeAnalytical
  end interface

  interface
     subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
       implicit none
       integer, intent(in) :: blkcnt
       integer, intent(in) :: blklst(blkcnt)
       integer, intent(in) :: nstep
       real, intent(in) :: dt
       real, intent(in) :: stime
     end subroutine Simulation_adjustEvolution
  end interface

  interface
     subroutine Simulation_freeUserArrays()
       implicit none

     end subroutine Simulation_freeUserArrays
  end interface

  interface
     subroutine Simulation_addStaticAcceleration(pos, sweepDir, blockID, numCells, grav)
       implicit none
       integer, dimension(2), intent(in) :: pos
       integer, intent(in)      :: sweepDir, blockID,  numCells
       real, intent(inout)  :: grav(numCells)
     end subroutine Simulation_addStaticAcceleration
  end interface

end Module Simulation_interface

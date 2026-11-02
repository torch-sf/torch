!!****h* source/Driver/Driver_interface
!!
!! NAME
!!  Driver_interface
!!
!! SYNOPSIS
!!
!!  use  Driver_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Driver module
!! that defines its public interfaces.
!!***
Module Driver_interface
 
  implicit none

# include "Flash.h"
# include "constants.h"

  interface
    subroutine Driver_abortFlash (errorMessage)
      implicit none
      character(len=*), intent(in) :: errorMessage
    end subroutine Driver_abortFlash
  end interface
  
  interface
    subroutine Driver_fatalError (errorMessage)
      implicit none
      character(len=*), intent(in) :: errorMessage
    end subroutine Driver_fatalError
  end interface
  
  interface
    subroutine Driver_dbgBreak()
    end subroutine Driver_dbgBreak
  end interface

  interface
    subroutine Driver_computeDt(nbegin, nstep, &
                        simTime, dtOld, dtNew)
      implicit none
      integer, intent(IN) :: nbegin, nstep
      real,    intent(IN) :: simTime    !! current simulation time
      real, intent(IN) :: dtOld      !! last time step we used
      real, intent(OUT):: dtNew      !! the new timestep we get. to be returned.
    end subroutine Driver_computeDt
  end interface

  interface
     subroutine Driver_computeCellLocations(nLocations, cellIsLocal, cellLocation, &
          nBlocks, blockList, doBroadcast)
       implicit none
       integer, intent(IN)  :: nLocations
       logical, intent(OUT) :: cellIsLocal(nLocations)
       integer, intent(OUT) :: cellLocation(5,nLocations)
       integer, intent(IN)  :: nBlocks
       integer, intent(IN)  :: blockList(nBlocks)
       logical, intent(IN),OPTIONAL  :: doBroadcast
     end subroutine Driver_computeCellLocations
  end interface

  interface
    subroutine Driver_driftBlock(src_file, src_line, blk, ptr, gds)
      implicit none
      character(len=*), intent(in) :: src_file
      integer, intent(in) :: src_line
      integer, intent(in) :: blk
      real, intent(in) :: ptr(:,:,:,:)
      integer, intent(in) :: gds
    end subroutine Driver_driftBlock
  end interface
  
  interface
    subroutine Driver_driftUnk(src_file, src_line, flags)
      implicit none
      character(len=*), intent(in) :: src_file
      integer, intent(in) :: src_line
      integer, intent(in), optional :: flags
    end subroutine Driver_driftUnk
  end interface
  
  interface
    subroutine Driver_evolveFlash()
      implicit none
    end subroutine Driver_evolveFlash
  end interface

  interface
    subroutine Driver_finalizeFlash()
      implicit none
    end subroutine Driver_finalizeFlash
  end interface

  interface
    subroutine Driver_finalizeSourceTerms( restart)
      implicit none
      logical, intent(in) :: restart
    end subroutine Driver_finalizeSourceTerms
  end interface

  interface
    subroutine Driver_getDt(dt)
      implicit none
      real, intent(out) :: dt
    end subroutine Driver_getDt
  end interface

  interface
    subroutine Driver_getElapsedWCTime(elapsedWCTime)
      implicit none
      real, intent(out) :: elapsedWCTime
    end subroutine Driver_getElapsedWCTime
  end interface

  interface
    subroutine Driver_getNStep(nstep)
      implicit none
      integer, intent(out) :: nstep
    end subroutine Driver_getNStep
  end interface

  interface
    subroutine Driver_getSimTime(simulationTime, simGeneration)
      implicit none
      real, intent(out) :: simulationTime
      integer, intent(out), OPTIONAL :: simGeneration
    end subroutine Driver_getSimTime
  end interface

  interface
    subroutine Driver_init()
      implicit none
    end subroutine Driver_init
  end interface

  interface
    subroutine Driver_initFlash()
      implicit none
    end subroutine Driver_initFlash
  end interface

  interface
    subroutine Driver_initMaterialProperties()
      implicit none
    end subroutine Driver_initMaterialProperties
  end interface

  interface
    subroutine Driver_initNumericalTools ()
    end subroutine Driver_initNumericalTools
  end interface

  interface
    subroutine Driver_initParallel ()
    end subroutine Driver_initParallel
  end interface

  interface
    subroutine Driver_initSourceTerms( restart)
      implicit none
      logical, intent(in) :: restart
    end subroutine Driver_initSourceTerms
  end interface

  interface
    subroutine Driver_sendOutputData()
      implicit none
    end subroutine Driver_sendOutputData
  end interface

  interface
    subroutine Driver_sourceTerms(blockCount, blockList, dt, pass)
      implicit none
      real, intent(IN)    :: dt
      integer, intent(IN) :: blockCount
      integer, dimension(blockCount), intent(IN):: blockList
      integer, OPTIONAL, intent(IN):: pass
    end subroutine Driver_sourceTerms
  end interface

  interface
    subroutine Driver_verifyInitDt()
      implicit none
    end subroutine Driver_verifyInitDt
  end interface

  interface
     subroutine Driver_getTimeStamp(dateStr)
       implicit none
       character(len=40), intent(OUT)     :: dateStr
     end subroutine Driver_getTimeStamp
  end interface

  interface
     subroutine Driver_putTimeStamp(dateStr)
       implicit none
       character(len=40), intent(IN)     :: dateStr
     end subroutine Driver_putTimeStamp
  end interface

  interface
     subroutine Driver_checkMPIErrorCode( errorCode)
       implicit none
       integer, intent(IN) :: errorCode
     end subroutine Driver_checkMPIErrorCode
  end interface


  interface
     subroutine Driver_superTimeStep(dt,nuSTS,nstepSTS,nstepTotalSTS,dt_subSTS)
       implicit none
       real, intent(IN)    :: dt,nuSTS
       integer, intent(IN) :: nstepSTS,nstepTotalSTS
       real, intent(OUT)   :: dt_subSTS
     end subroutine Driver_superTimeStep
  end interface


  interface
     subroutine Driver_getMype(communicatorType, mype, axis)
       integer, INTENT(IN) :: communicatorType
       integer, INTENT(OUT) :: mype
       integer, optional, intent(IN) :: axis
     end subroutine Driver_getMype
  end interface

  interface
     subroutine Driver_getNumProcs(communicatorType, numProcs, axis)
       integer, INTENT(IN) :: communicatorType
       integer, INTENT(OUT) :: numProcs
       integer, optional, intent(IN) :: axis
     end subroutine Driver_getNumProcs
  end interface

  interface
     subroutine Driver_getComm(communicatorType, communicator, axis)
       integer, INTENT(IN) :: communicatorType
       integer, INTENT(OUT) :: communicator
       integer, optional, intent(IN) :: axis
     end subroutine Driver_getComm
  end interface

  interface
     subroutine Driver_setupParallelEnv ()
       
       implicit none             
     end subroutine Driver_setupParallelEnv
  end interface

  interface
     subroutine Driver_mpiThreadSupport(mpiThreadSupport)
       implicit none
       logical, intent(OUT) :: mpiThreadSupport
     end subroutine Driver_mpiThreadSupport
  end interface

  interface
     subroutine Driver_logMemoryUsage(callsite)
       implicit none
       character(len=*),intent(IN) :: callsite
     end subroutine Driver_logMemoryUsage
  end interface
  
  interface 
     subroutine Driver_diagnostics (blockCount, blockList, dt)
        implicit none
        real,    intent(IN) :: dt
        integer, intent(IN) :: blockCount
        integer, dimension(blockCount), intent(IN):: blockList
     end subroutine Driver_diagnostics
  end interface 

end Module Driver_interface

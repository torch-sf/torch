!!
!! NAME
!!  
!!  rt_mean_intensity_block
!!
!!
!! SYNOPSIS
!! 
!!  call rt_mean_intensity_block(block_no)
!!  call rt_mean_intensity_block(integer)
!!  
!! DESCRIPTION
!!
!!  
!!
!!***
subroutine rt_mean_intensity_block (block_no)
  !
#include "Flash.h"
#include "constants.h"
  
  !
  use Grid_interface, ONLY             : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY           : Driver_abortFlash      
  use RadTrans_hybridCharModule,  ONLY : dMeanMax
  use rt_data, ONLY                    : rt_ALI, rt_epsilon
  !
  implicit none
  !
  include 'Flash_mpi.h'
  !
  integer, INTENT(in)               :: block_no
  !
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  !
  integer :: blkLimitsGC(LOW:HIGH,MDIM), blkLimits(LOW:HIGH,MDIM)
  integer :: i, j, k
  ! 
  real :: meanFS,meanOld,meanNew,dMean 
  real :: lambda,lambdaInv 
  !

  call Grid_getBlkPtr(block_no, solnData)
  call Grid_getBlkIndexLimits(block_no,blkLimits,blkLimitsGC)
      
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           
           !
           ! get the formal solution
           !
           meanFS       = solnData(MEAN_VAR,i,j,k)
           meanOld      = solnData(MEAO_VAR,i,j,k)
           !
           !  now comes the (accelerated) lambda iteration step
           !
           select case(rt_ALI)
              !
           case(0)
              !
              ! this is the ordinary lambda iteration without acceleration
              ! it is trivial, the new radiation field is simply the
              ! formal solution.
              !
              meanNew    = meanFS
              !
              dMean = abs(meanNew-meanOld)/MAX(meanOld,1.e-100)
           case(1)
              !
              ! this is the lambda iteration with a local lambda star operator.
              ! For rt_epsilon = 1. this reduces to case(0), e.g. ordinary lambda
              ! iteration.
              !
              lambda     = solnData(LAMB_VAR,i,j,k)
              lambdaInv  = 1.d0/(1.d0-lambda*(1.-rt_epsilon))
              IF(meanOld.NE.0.) THEN
                meanNew = lambdaInv * (meanFS - meanOld)
                dMean = abs(meanNew)/MAX(meanOld,1.e-100)
                meanNew = meanNew + meanOld
              ELSE
                meanNew = meanFS
                dMean = abs(meanNew-meanOld)/MAX(meanOld,1.e-100)
              END IF
              !
           end select

           if(meanNew < 0.0) then
              write(*,*) '========================================='
              write(*,*) 'ERROR in the Lambda Iteration'
              write(*,*) 'meanNew is negative!'
              write(*,*) 'meanNew:', meanNew
              write(*,*) 'meanFS:', meanFS
              write(*,*) 'meanOld:', meanOld
              write(*,*) 'lambda:', lambda
              write(*,*) 'i,j,k:', i,j,k
              !write(*,*) 'aborting flash...'
              write(*,*) 'setting meanNew=meanFS'
              write(*,*) '========================================='
              call Driver_abortFlash('radiation: ERROR raytrace_3DRT')
           endif
           !
           if(dMean>dMeanMax) dMeanMax=dMean
           !
           solnData(MEAN_VAR,i,j,k) = meanNew
           !
        enddo
     enddo
  enddo
  !
  call Grid_releaseBlkPtr(block_no, solnData)
  !
  return
  !
end subroutine rt_mean_intensity_block
!

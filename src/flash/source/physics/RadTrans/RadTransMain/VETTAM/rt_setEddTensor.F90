!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/rt_setEddTensor
!!
!! NAME
!!  
!!  rt_setEddTensor
!!
!!
!! SYNOPSIS
!! 
!!  rt_setEddTensor()
!!  
!! DESCRIPTION
!!
!!	Set the Eddington tensor
!!
!!***
!!***

#include "Flash.h"
#include "constants.h"
#include "Multispecies.h"
SUBROUTINE rt_setEddTensor()
  use RadTrans_data
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_getBlkIndexLimits, Grid_releaseBlkPtr, Grid_getCellCoords
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
#ifdef RAYTRACE_3DRT
  use raytrace_data, ONLY: rt_useRayTrace
#endif
  real :: fred, chi
  integer :: sizeX, sizeY, sizeZ, istat
  real, dimension(:), allocatable :: x, y, z
  integer :: b, blockID, myPE, i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData
  real :: xhat,yhat,zhat,mohmag, r
  logical :: gcell=.true.

  SELECT CASE(rt_etens)

  !Eddington Approximation
  CASE('edd','fld')

    do b = 1, blockCount
      
      blockID = blockList(b)
      call Grid_getBlkPtr(blockID,solnData)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      solnData(XXED_VAR,:,:,:) = 1./3.
#if NDIM>1
      solnData(YYED_VAR,:,:,:) = 1./3.
      solnData(XYED_VAR,:,:,:) = 0.0
#if NDIM>2
      solnData(ZZED_VAR,:,:,:) = 1./3.
      solnData(YZED_VAR,:,:,:) = 0.0
      solnData(XZED_VAR,:,:,:) = 0.0
#endif
#endif
      call Grid_releaseBlkPtr(blockID,solnData)
    end do

  !Moment 1 condition (time-lagged)
  CASE('m1','mom1')

    do b = 1, blockCount
      
      blockID = blockList(b)
      call Grid_getBlkPtr(blockID,solnData)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
        do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
          do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
            ! Calculate reduced flux
            mohmag = solnData(MOHX_VAR,i,j,k)**2 
#if NDIM>1
            mohmag = mohmag + solnData(MOHY_VAR,i,j,k)**2 
#if NDIM>2
            mohmag = mohmag + solnData(MOHZ_VAR,i,j,k)**2 
#endif
#endif
            mohmag = SQRT(mohmag)
            fred = mohmag/(rt_speedlt*solnData(ERAD_VAR,i,j,k))
            if(solnData(ERAD_VAR,i,j,k) .eq. 0.0) fred = 0.0

            !Safety check
            if(fred .gt. huge(fred)) then
              call Driver_abortFlash("Reduced flux is infinite.")
            endif

            if(fred .gt. 1) then 
              fred = 1.0
            endif
            
            !Calculate Eddington Factor from Levermore (1984) relation
            chi = (3+4*fred**2)/(5+2*SQRT(4.-3.*fred**2))

            if((1./mohmag) .gt. huge(mohmag)) then 
              xhat = 0.
              yhat = 0.
              zhat = 0.
            else
              !Calculate unit vectors
              xhat = solnData(MOHX_VAR,i,j,k)/mohmag
              yhat = solnData(MOHY_VAR,i,j,k)/mohmag
              zhat = solnData(MOHZ_VAR,i,j,k)/mohmag
            endif

            !Now set Tensor based on Levermore conditions
            solnData(XXED_VAR,i,j,k) = (1-chi)/2. + (3*chi-1)/2. * xhat * xhat
#if NDIM>1
            solnData(YYED_VAR,i,j,k) = (1-chi)/2. + (3*chi-1)/2. * yhat * yhat
            solnData(XYED_VAR,i,j,k) = (3*chi-1)/2. * xhat * yhat
#if NDIM>2
            solnData(ZZED_VAR,i,j,k) = (1-chi)/2. + (3*chi-1)/2. * zhat * zhat
            solnData(YZED_VAR,i,j,k) = (3*chi-1)/2. * yhat * zhat
            solnData(XZED_VAR,i,j,k) = (3*chi-1)/2. * xhat * zhat
#endif
#endif

          end do
        end do
      end do
      call Grid_releaseBlkPtr(blockID,solnData)
    end do

  ! Calculate Etensors with the ray-tracer
  CASE('vet')

#ifdef RAYTRACE_3DRT
    ! Call Ray-tracer to compute instantaneous eddington factors
    if(.not. rt_useRayTrace) call Driver_abortFlash("rt_etens = vet requires useRayTrace = true.")
    call Timers_start("vet_raytrace")
    call Logfile_stamp("Computing Eddington Tensor", "[VET] ")
    call RadTrans_RayTrace_Diffuse()
    call Timers_stop("vet_raytrace")
#else
    call Driver_abortFlash("rt_etens = vet requires the 3DRT to be compiled in.")
#endif
  

  ! Set a radial eddington tensor (appropriate for an optically thin medium with a central source)
  ! This is effectively a tensor with the r^r^ term non-zero, and the rest zero, translated to Cartesian coordinates
  CASE('radial')
    do b = 1, blockCount
      blockID = blockList(b)
      call Grid_getBlkPtr(blockID,solnData)
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
      sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
      sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
      sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1
      allocate(x(sizeX), stat=istat)
      allocate(y(sizeY), stat=istat)
      allocate(z(sizeZ), stat=istat)
      call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, z, sizeZ)
      call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, y, sizeY)
      call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, x, sizeX)

      do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
        do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
          do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
            r = (x(i))**2
#if NDIM>1
            r = r+ (y(j))**2
#if NDIM>2
            r = r + (z(k))**2
#endif
#endif
            r = SQRT(r)
            solnData(XXED_VAR,i,j,k) = (x(i)**2)/(r**2)
#if NDIM>1
            solnData(YYED_VAR,i,j,k) = (y(j)**2)/(r**2)
            solnData(XYED_VAR,i,j,k) = (x(i)*y(j))/(r**2)
#if NDIM>2
            solnData(YZED_VAR,i,j,k) = (y(j)*z(k))/(r**2)
            solnData(XZED_VAR,i,j,k) = (x(i)*z(k))/(r**2)
#endif
#endif

          end do
        end do
      end do
      call Grid_releaseBlkPtr(blockID,solnData)
    end do
    deallocate(x)
    deallocate(y)
    deallocate(z)

  ! Here the user can implement their own version, or set something to a fixed value in the initial conditions
  CASE('user')
    !Nothing to be done here

  CASE DEFAULT
    call Driver_abortFlash("[VET] : Unrecognised eddington tensor type.") 
  END SELECT

end SUBROUTINE rt_setEddTensor
!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/rt_sinkHydro
!!
!! NAME
!!  
!!  rt_sinkHydro
!!
!!
!! SYNOPSIS
!! 
!!  rt_sinkHydro()
!!  
!! DESCRIPTION
!!
!!	Set the Eddington tensor
!!
!!***
!!***
SUBROUTINE rt_sinkHydro(blockCount_,blockList_,dt)

  use RadTrans_data
  use Hydro_data, ONLY: hy_useHydro
  use Grid_interface, ONLY: Grid_getBlkIndexLimits,Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface, ONLY:Eos_wrapped
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, INTENT(IN) :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_  
  real, intent(in) :: dt
  real             :: rho,temp,vx,vy,vz
  real             :: dE,dMomx,dMomy,dMomz, ekin_old, ekin_new, d_ekin
  real             :: cellsize(MDIM)
  integer          :: b,blockID
  integer          :: i,j,k,l 
  real             :: d1x,d1y,d1z,d2,opac_planck,opac_rosseland 
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,cellsize)
    do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
      do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

          rho    = solnData(DENS_VAR,i,j,k)
#ifdef TEMP_VAR
          temp   = solnData(TEMP_VAR,i,j,k)
#endif
#ifdef VELX_VAR 
          vx = solnData(VELX_VAR,i,j,k)
#endif          
#ifdef VELY_VAR
          vy = solnData(VELY_VAR,i,j,k)
#endif
#ifdef VELZ_VAR
          vz = solnData(VELZ_VAR,i,j,k)
#endif          
          d1x = dt/(2*cellsize(IAXIS))
          d1y = dt/(2*cellsize(JAXIS))
          d1z = dt/(2*cellsize(KAXIS))
          d2 = dt * rt_speedlt
          opac_planck = solnData(TAUP_VAR,i,j,k) 
          opac_rosseland = solnData(TAUR_VAR,i,j,k) 

          !Reset to zero
          dMomx = 0.0
          dMomy = 0.0
          dMomz = 0.0
          ekin_old = 0.0
          ekin_new = 0.0
          dE = 0.0

          ! Add Sink Heating
#ifdef FSPT_VAR
          if(.not. rt_sink_implicit .and. rt_sinkheat) then
            dE = dE + solnData(FSPT_VAR,i,j,k)*dt
          else
            dE = 0.0
          endif
#endif

          ekin_old = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
          ekin_old = ekin_old+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
          ekin_old = ekin_old + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif

          ! Add momentum only if hydro switched on
          if(hy_useHydro .and. rt_sinkmom) then

          !Radiation Pressure from stellar radiation if any
#ifdef STMX_VAR
            dMomx = dMomx + solnData(STMX_VAR,i,j,k)*dt
#if NDIM>1
            dMomy = dMomy + solnData(STMY_VAR,i,j,k)*dt
#if NDIM>2        
            dMomz = dMomz + solnData(STMZ_VAR,i,j,k)*dt
#endif 
#endif        

#endif        
          endif

#ifdef EINT_VAR
          solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + dE/solnData(DENS_VAR,i,j,k)
#else
          solnData(ERAD_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k) + dE
#endif

          ! Update Gas velocities now. 

          if(solnData(DENS_VAR,i,j,k) .gt. tiny(solnData(DENS_VAR,i,j,k))) then
        
            solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + dMomx/solnData(DENS_VAR,i,j,k)
#if NDIM>1
            solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + dMomy/solnData(DENS_VAR,i,j,k)
#if NDIM>2
            solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + dMomz/solnData(DENS_VAR,i,j,k)  
#endif 
#endif          
          else 
            solnData(VELX_VAR,i,j,k) = 0.0
            solnData(VELY_VAR,i,j,k) = 0.0
            solnData(VELZ_VAR,i,j,k) = 0.0
          endif
        
          ekin_new = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
          ekin_new = ekin_new+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
          ekin_new = ekin_new + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif

          ! compute energy injection
          d_ekin = ekin_new - ekin_old

#ifdef ENER_VAR
          ! update the total energy to be consistent with new velocities
          if(solnData(DENS_VAR,i,j,k) .gt. tiny(solnData(DENS_VAR,i,j,k))) then
            solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + d_ekin + dE/solnData(DENS_VAR,i,j,k)
          else
            solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)
          endif
#endif
       
        end do
      end do
    end do 
    call Grid_releaseBlkPtr(blockID,solnData)
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
  end do   

END SUBROUTINE rt_sinkHydro



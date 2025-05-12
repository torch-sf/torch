!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/rt_dustTerms
!!
!! NAME
!!  
!!  rt_dustTerms
!!
!!
!! SYNOPSIS
!! 
!!  rt_dustTerms()
!!  
!! DESCRIPTION
!!
!!	Imparts the dust radiation pressure and any reprocessing of UV to IR wavelengths
!!
!!***
!!***

SUBROUTINE rt_dustTerms(dt)
  use RadTrans_data
  use Grid_interface, ONLY: Grid_getBlkIndexLimits,Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface, ONLY:Eos_wrapped
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none


#include "Flash.h"
#include "constants.h"

    
  real, intent(in) :: dt
  real             :: rho,temp,vx,vy,vz
  real             :: dE,dMomx,dMomy,dMomz, ekin_old, ekin_new, d_ekin
  real             :: cellsize(MDIM)
  integer          :: b,blockID
  integer          :: i,j,k,l 
  real             :: d1x,d1y,d1z,d2,opac_planck,opac_rosseland 
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData

  do b = 1, blockCount
    blockID = blockList(b)
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
          !Planck and Rosseland opacity of dust; this is obtained by subtracting out the opacity due to hydrogen
          opac_planck = solnData(TAUP_VAR,i,j,k)
          opac_rosseland = solnData(TAUR_VAR,i,j,k)
          !Subtract hydrogen contribution for EUV band (since this is added by the photoionization module)
#ifdef TAUH_VAR
          if(current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
          & current_band .eq. 'EUV_15P2_INFTY') then
            opac_planck = opac_planck - solnData(TAUH_VAR,i,j,k)
            opac_rosseland = opac_rosseland - solnData(TAUH_VAR,i,j,k)
          endif
#endif

#ifdef REIR_VAR
          !Accumulate the energy absorbed by dust in all bands by assuming it is re-emitted in the IR completely
          solnData(REIR_VAR,i,j,k) = solnData(REIR_VAR,i,j,k) + opac_planck * rt_speedlt * solnData(ERAD_VAR,i,j,k) * dt
#endif
          !Photoelectric band: the absorbed energy is stored since it will be used by KROME for PE heating of gas
#ifdef PEFL_VAR
          if(current_band .eq. 'FUV') then 
          !Accumulate the energy absorbed by dust in the PE band
            solnData(PEFL_VAR,i,j,k) = solnData(PEFL_VAR,i,j,k) + opac_planck * rt_speedlt * solnData(ERAD_VAR,i,j,k) * dt &
                    & * rt_speedlt/1.6e-3 
          endif
#endif

          !LW band: the absorbed energy (by dust) is stored since it will be used for photodiassociation of H2, CO and C ionization; and for PE heating by dust
#ifdef ULWD_VAR
          if(current_band .eq. 'LW' .or. current_band .eq. 'LYMAN_WERNER') then 
          !Accumulate the energy absorbed by dust in the LW band by assuming it is re-emitted in the IR completely
            solnData(ULWD_VAR,i,j,k) = solnData(ULWD_VAR,i,j,k) + opac_planck * rt_speedlt * solnData(ERAD_VAR,i,j,k) * dt
          endif
#endif

          !Reset to zero
          dMomx = 0.0
          dMomy = 0.0
          dMomz = 0.0
          ekin_old = 0.0
          ekin_new = 0.0

          ekin_old = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
          ekin_old = ekin_old+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
          ekin_old = ekin_old + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif
                
          !Computing explicit terms to add to the gas momentum

          !Gas-Radiation momentum exchange term
          dMomx = dMomx + dt * opac_rosseland * solnData(MOHX_VAR,i,j,k)/rt_speedlt
#if NDIM>1
          dMomy = dMomy + dt * opac_rosseland * solnData(MOHY_VAR,i,j,k)/rt_speedlt
#if NDIM>2
          dMomz = dMomz + dt * opac_rosseland * solnData(MOHZ_VAR,i,j,k)/rt_speedlt
#endif
#endif

          !Only add these terms if VET switched on
          if(rt_ovcterms) then
          
            dMomx = dMomx - &
              d2*opac_rosseland *vx*(1+solnData(XXED_VAR,i,j,k))*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
#if NDIM>1

            dMomx = dMomx - d2*opac_rosseland*vy*solnData(XYED_VAR,i,j,k)*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2 

            dMomy = dMomy - &
              d2*opac_rosseland * (vy*(1+solnData(YYED_VAR,i,j,k)) + &
                vx*solnData(XYED_VAR,i,j,k)) * solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
#if NDIM>2


            dMomx = dMomx - d2*opac_rosseland*vz*solnData(XZED_VAR,i,j,k)*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
            dMomy = dMomy - d2*opac_rosseland*vz*solnData(YZED_VAR,i,j,k)*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2

            dMomz = dMomz - &
              d2*opac_rosseland * (vz*(1+solnData(ZZED_VAR,i,j,k)) + vx*solnData(XZED_VAR,i,j,k) + &
              vy*solnData(YZED_VAR,i,j,k)) * solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
#endif
#endif
          endif

        ! Update Gas velocities now. 
        
        solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + dMomx/solnData(DENS_VAR,i,j,k)
#if NDIM>1
        solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + dMomy/solnData(DENS_VAR,i,j,k)
#if NDIM>2
        solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + dMomz/solnData(DENS_VAR,i,j,k)  
#endif 
#endif

        
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
        solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + d_ekin
#endif
       
        end do
      end do
    end do 
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do 


END SUBROUTINE rt_dustTerms

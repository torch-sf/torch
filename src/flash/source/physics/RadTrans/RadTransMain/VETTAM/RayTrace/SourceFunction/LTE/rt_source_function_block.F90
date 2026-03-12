!!
!! NAME
!!  
!!  source_function_block
!!
!!
!! SYNOPSIS
!! 
!!  call source_function_block( block_no)
!!  call source_function_block(integer)
!!  
!! DESCRIPTION
!!
!!  
!!
!!***

subroutine rt_source_function_block (block_no)

#include "Flash.h"
#include "constants.h"
        
      use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
      use RadTrans_hybridCharModule,  ONLY: dSourceMax, sb_law
      use raytrace_data, ONLY: rt_epsilon
      use RadTrans_data, ONLY: current_band, rt_radconst
#ifdef IONY_MSCALAR
      use Eos_data, only: eos_singleSpeciesA
      use PhysicalConstants_interface, ONLY: PhysicalConstants_get
      use rt_ionisedata, ONLY: alpha_ground_constant, ion_ots, alpha_type
      use rt_ionisemodule, ONLY: get_recombination_ground
#elif defined(IHA_SPEC) && defined(UEUV_VAR)
      use rt_ionisedata, ONLY : hpA, elecA, alpha_ground_constant, ion_ots, alpha_type
      use rt_ionisemodule, ONLY: get_recombination_ground
#endif

      implicit none

#include "Flash_mpi.h"

      real, DIMENSION(:,:,:,:), POINTER :: solnData
      integer, INTENT(in) :: block_no

      integer :: blkLimitsGC(LOW:HIGH,MDIM), blkLimits(LOW:HIGH,MDIM)
      integer :: i, j, k

      real :: tmp
      real :: source,sourceOld,mean
      real :: dSource 
      real :: t0, cv_inv, jstr, taur, diffuse_euv_j
      real :: mH, iony, nH, hnu, ne , nHplus, temp_elec, alpha_ground_cell

      call Grid_getBlkPtr(block_no, solnData)
      call Grid_getBlkIndexLimits(block_no,blkLimits,blkLimitsGC)
      
      do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
         do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
            do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              !Use gas temperature if available. If dust temperature stored in a different variable don't use it.
#if defined(TEMP_VAR) && !defined(KDUS_VAR)
               tmp       = solnData(TEMP_VAR,i,j,k)
#else
               !In case of Isothermal-like EOS or condition of radiative eq, use radiation temperature
               tmp       = (solnData(ERAO_VAR,i,j,k)/rt_radconst)**0.25
#endif
               sourceOld = solnData(SOUR_VAR,i,j,k)
               mean      = solnData(MEAN_VAR,i,j,k)
               if(current_band .eq. 'IR') source = (1. - rt_epsilon) * mean + rt_epsilon * sb_law(tmp)

#ifdef ERAD_VAR
               !The emission term (from gas) is zero for the UV bands
               if(current_band .eq. 'FUV' .or. current_band .eq. 'EUV' .or. &
               current_band .eq. 'LYMAN_WERNER' .or. current_band .eq. 'LW' .or. &
               current_band .eq. 'PE' .or. current_band .eq. 'EUV_13P6_15PE' .or. &
               current_band .eq. 'EUV_15P2_INFTY') source = 0.0
#endif

               ! Add contribution from diffuse emission due to recombinations to ground state
#ifdef UEUV_VAR
               if(current_band .eq. 'EUV' .and. .not. ion_ots) then 
                  source = 0.0
                  taur = solnData(TAUR_VAR,i,j,k)
                  !Hardcoded now; change!
                  hnu      = 13.6*1.6021764620000066e-12
                  !Add diffuse EUV recombination field. alpha_ground will be >0 if ion_ots is .false.
                  !Get the coefficient for recombinations to ground state
#ifdef TEMP_VAR
                  temp_elec = solnData(TEMP_VAR,i,j,k)
#elif defined(TGAS_VAR)
                  temp_elec = solnData(TGAS_VAR,i,j,k)
#endif
                  if(alpha_type .eq. 'default') then
                        alpha_ground_cell = get_recombination_ground(temp_elec,ion_ots)
                  else if(alpha_type .eq. 'constant') then
                        alpha_ground_cell = alpha_ground_constant
                  endif
#ifdef IONY_MSCALAR
                  call PhysicalConstants_get("proton mass",mH)
                  iony = solnData(IONY_MSCALAR,i,j,k)
                  nH = solnData(DENS_VAR,i,j,k)/(eos_singleSpeciesA*mH)
                  nHplus = nH*(1.-iony)
                  ne = nHplus

                  diffuse_euv_j = alpha_ground_cell * ne * nHplus * hnu/(4*PI*taur)
                  
                  !KROME Version
#elif defined(IHA_SPEC) && defined(UEUV_VAR)
                  nHplus = solnData(DENS_VAR,i,j,k)*solnData(IHP_SPEC,i,j,k)/hpA
                  ne     = solnData(DENS_VAR,i,j,k)*solnData(ELEC_SPEC,i,j,k)/(elecA)
                  diffuse_euv_j = alpha_ground_cell * ne * nHplus * hnu/(4*PI*taur)
#else
                  diffuse_euv_j = 0.0 !No source of EUV photons if these variables not defined
                  
#endif
                  !Note/TODO: since the RT uses the source function definition, which divides by specific opacity, it can cause problems if tau=0
                  !This can be a problem for fully ionised material (iony=1, and therefore taur=0); one fix might be to set a floor on the opacity.
                  !TODO: DO THIS AT SOME POINT; LEAVING IT HERE FOR NOW.
                  if(taur .gt. 0) &
                  & source = source + diffuse_euv_j
                  
               endif
#endif

#ifdef JSTR_VAR
               jstr = solnData(JSTR_VAR,i,j,k)
               taur = solnData(TAUR_VAR,i,j,k)
               if(taur .gt. 0) &
               & source = source + jstr/(4*PI*taur)
#endif

               dSource = abs(source-sourceOld)/(sourceOld+1.d99)
               !if (dSource.gt.dSourceMax) dSourceMax = dSource

               solnData(SOUR_VAR,i,j,k) = source 

            enddo
         enddo
      enddo

      call Grid_releaseBlkPtr(block_no, solnData)
      return

end subroutine rt_source_function_block

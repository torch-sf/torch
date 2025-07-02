
!!****f* source/physics/RadTrans/RadTransMain/VETTAM/Photoionization/OnlyH/rt_ionRateEqs
!!
!! NAME
!!  
!!  rt_ionRateEqs
!!
!!
!! SYNOPSIS
!! 
!!  rt_ionRateEqs(integer(IN) :: blockCount
!!       integer(IN) :: blockList(blockCount),
!!          real(IN) :: dt)
!!  
!! DESCRIPTION
!!
!! Solves the rate-equation for the ionisation state of atomic H gas provided
!! the EUV radiation field and number densities of species are known.
!! Currently uses a simple, explicit, analytical prescription for the 
!! ionisation state. Future versions may include implicit schemes
!! and/or other species. 
!!
!!
!! AUTHOR
!!  Shyam Harimohan Menon (2022-2023)
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!! 
!!
!!***
#include "Flash.h"
#include "constants.h"
#include "Multispecies.h"
#ifdef H2_SPEC
SUBROUTINE IonizeH2(blockCount_, blockList_, dt, time)
  use rt_ionisedata
  use Grid_interface
  use Grid_data, ONLY: gr_smallx
  use RadTrans_data, ONLY: current_band
  implicit none
  integer, intent(IN)                         :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                            :: dt,time
  integer                                     :: b,blockID, i, j, k, n
  integer,dimension(2,MDIM)                   :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer             :: solnData
  real                                        :: nH0, nH2, sum, frac_15p2_infty_H2, frac_15p2_infty_H, nH2_old, dnH2
  real                                        :: IonizationRate_H2, IonizationRate_H, IonizationRate, Gamma_H2

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)


          nH0 = solnData(DENS_VAR,i,j,k)*solnData(IHA_SPEC,i,j,k)/hA
          nH2 = solnData(DENS_VAR,i,j,k) * solnData(H2_SPEC,i,j,k)/h2A
          IonizationRate    = solnData(UEUV_VAR,i,j,k)/(hnu)
          !If H2 ionization is switched on, compute the fraction of photons absorbed in the hnu>15.2 eV band by H and H2 respectively
          if(current_band .eq. 'EUV_15P2_INFTY' .and. useH2Ionize) then
            frac_15p2_infty_H2 = ion_sigmaH2_15p2_infty*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH_15p2_infty*nH0/(solnData(TAUH_VAR,i,j,k))
          !If there is only one EUV band >13.6 eV
          else if(current_band .eq. 'EUV' .and. useH2Ionize) then
            frac_15p2_infty_H2 = ion_sigmaH2*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH*nH0/(solnData(TAUH_VAR,i,j,k))
          else if(current_band .eq. 'EUV_13P6_15P2' .or. .not. useH2Ionize) then
            frac_15p2_infty_H2 = 0.0
            frac_15p2_infty_H = 1.0
          endif

          IonizationRate_H2 = IonizationRate * frac_15p2_infty_H2
          IonizationRate_H = IonizationRate * frac_15p2_infty_H

          if(useH2Ionize) then
            nH2_old = nH2
            !dnH2 = IonizationRate_H2*dt
            Gamma_H2 = IonizationRate_H2/nH2
            if(nH2 .eq. 0) Gamma_H2 = 0.0
            dnH2 = (1.0 - exp(-Gamma_H2 * dt))*nH2_old
            nH2 = nH2_old - dnH2
            if(nH2 .lt. -1.e-5) then
              print *, 'nH2 is less than zero. nH2,nH2old, dnH2 ', nH2, nH2_old, dnH2
            else
              nH2 = max(nH2,0.0)
            endif
            !Convert to mass fraction
            solnData(H2_SPEC,i,j,k) = nH2 * h2A/solnData(DENS_VAR,i,j,k)

            !Now add contributions from ionized molecular H2
            nH0 = nH0 + 2*dnH2
            solnData(IHA_SPEC,i,j,k) = nH0 * hA/solnData(DENS_VAR,i,j,k)

            !Renormalise
            sum = 0.0
            do n = SPECIES_BEGIN, SPECIES_END
              sum = sum + max(gr_smallx, solnData(n,i,j,k))
            enddo

            ! re-normalise sum of species fractions to 1
            do n = SPECIES_BEGIN, SPECIES_END
              solnData(n,i,j,k) = solnData(n,i,j,k) / sum
            enddo

          endif
          
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

END SUBROUTINE IonizeH2
#endif

SUBROUTINE IonizeRateEquation(blockCount_,blockList_,dt,time)
  use rt_ionisedata
  use rt_ionisemodule, ONLY: get_recombination_coefficient
  use Grid_interface
  use Driver_data, ONLY: dr_globalMe
  use Eos_data, ONLY: eos_singleSpeciesA
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Grid_data, ONLY: gr_smallx
  use RadTrans_data, ONLY: current_band
  implicit none

! #include "Flash.h"
! #include "constants.h"
! #include "Multispecies.h"

  integer, intent(IN)                         :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                            :: dt,time
  integer                                     :: b,blockID, i, j, k, n
  integer,dimension(2,MDIM)                   :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer             :: solnData
  real                                        :: rho, iony, iono, nH, nH0, nHplus, ne, ne_new, sum, temp, alpha_rec_cell
  real, save                                  :: mH, pi, kb
  logical, save                               :: first_call = .true.
  real                                        :: RecombinationRate, IonizationRate, Gamma, xeq_0, x0, tir, xn, xplus, Kvar, Gamma_H2
  real                                        :: IonizationRate_H, IonizationRate_H2, dnH_by_dt, nH2
  real                                        :: frac_15p2_infty_H2, frac_15p2_infty_H, nH2_old, dnH2

  if(first_call) then
    call PhysicalConstants_get("proton mass",mH)
    call PhysicalConstants_get("pi",pi)
    call PhysicalConstants_get("Boltzmann",kb)
    first_call = .false.
  endif
  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)

    do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
      do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i = blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)

          rho    = solnData(DENS_VAR,i,j,k)
#ifdef IHA_SPEC
          ! Neutral fraction
          !Neutral fraction = nH0/(nH0 + nHplus)
          ! Total hydrogen nuclei number density = nHplus + nH0
          nH0 = rho*solnData(IHA_SPEC,i,j,k)/hA
          nHplus = rho*solnData(IHP_SPEC,i,j,k)/hpA
          nH     = nH0 + nHplus
          ! ne = mass fraction * density/m_e
          ne       = solnData(ELEC_SPEC,i,j,k) * rho/(elecA) 
          !Neutral fraction = nH0/(nH0 + nHplus)
          iony = nH0/nH
#elif defined(IONY_MSCALAR)
          nH = rho/(eos_singleSpeciesA*mH)
          iony   = solnData(IONY_MSCALAR,i,j,k)
          nH0 = iony * nH
          nHplus = nH * (1.-iony)
          nHplus = MAX(MIN(nHplus,nH),0.0) !Prevent unphysical values for nHplus
          ! Electron number density
          ne     = nHplus
#endif
          iony   = MAX(MIN(iony,1.0),0.0)
          ! Old Neutral fraction
          iono   = solnData(IONO_VAR,i,j,k)
          ! Number density of neutral hydrogen
#ifdef UEUV_VAR
          !Volumetric ionization rate
          IonizationRate    = solnData(UEUV_VAR,i,j,k)/(hnu)
#else
          !NOTE: If ray-tracer not included then user can define something here
          IonizationRate     = 0.0
#endif

#if defined(H2_SPEC)
          nH2 = rho * solnData(H2_SPEC,i,j,k)/h2A
          !If H2 ionization is switched on, compute the fraction of photons absorbed in the hnu>15.2 eV band by H and H2 respectively
          if(current_band .eq. 'EUV_15P2_INFTY' .and. useH2Ionize) then
            frac_15p2_infty_H2 = ion_sigmaH2_15p2_infty*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH_15p2_infty*nH0/(solnData(TAUH_VAR,i,j,k))
          !If there is only one EUV band >13.6 eV
          else if(current_band .eq. 'EUV' .and. useH2Ionize) then
            frac_15p2_infty_H2 = ion_sigmaH2*nH2/(solnData(TAUH_VAR,i,j,k))
            frac_15p2_infty_H = ion_sigmaH*nH0/(solnData(TAUH_VAR,i,j,k))
          else if(current_band .eq. 'EUV_13P6_15P2' .or. .not. useH2Ionize) then
            frac_15p2_infty_H2 = 0.0
            frac_15p2_infty_H = 1.0
          endif

          IonizationRate_H2 = IonizationRate * frac_15p2_infty_H2
          IonizationRate_H = IonizationRate * frac_15p2_infty_H
            
#else       
          !Pure Hydrogen case
          IonizationRate_H2 = 0.0
          IonizationRate_H = IonizationRate
#endif
          !Get recombination coefficient
#ifdef TEMP_VAR
          temp = solnData(TEMP_VAR,i,j,k)
#elif defined(TGAS_VAR)
          temp = solnData(TGAS_VAR,i,j,k)
#endif
          if(alpha_type .eq. 'default') then
            alpha_rec_cell = get_recombination_coefficient(temp,ion_ots)
          else if(alpha_type .eq. 'constant') then
            alpha_rec_cell = alpha_rec_constant
          else
            print *, alpha_type
            call Driver_abortFlash("[RateEquations.F90]:alpha_type should be 'constant' or 'default'; check!")
          endif

          if(multiple_ionbands) then
            IonizationRate_H = IonizationRate_H + solnData(HEUV_VAR,i,j,k)/(hnu)
          endif

          !Recombination rates
          RecombinationRate = alpha_rec_cell * ne * nHplus
          !Ionization rate
          Gamma = IonizationRate_H/(nH0)

          !Prevent NaN in Gamma if completely ionized cell
          if(nH0 .eq. 0) &
            & Gamma = 0.0

#ifdef H2_SPEC
          if(ion_type .ne. 1 .and. useH2Ionize) call Driver_abortFlash("[VET_Ionize]: H2 ionization is only supported with &
          & ion_type =1. Change ion_type or implement yourself :P")
#endif
          
          !1: Analytical explicit solution assuming \alpha_B, Gamma, and n_e constant in update
          if(ion_type .eq. 1) then

            xeq_0 = (alpha_rec_cell*ne)/(Gamma+alpha_rec_cell*ne)
            x0    = iony
            tir   = (Gamma+alpha_rec_cell*ne)**(-1)
            if(IonizationRate_H .eq. 0) then 
              !If fully neutral, remain neutral
              if(ne .eq. 0) then
                xn = x0 
              else
                !Eq 19 of Kim et al 2017
                xn = (x0+(1-x0)*alpha_rec_cell*nH*(time+dt))/(1+(1-x0)*alpha_rec_cell*nH*(time+dt))
              endif

            else
              ! Eq 20 of Kim et al 2017
              xn    = xeq_0 + (x0 - xeq_0)*exp(-dt/tir)
            endif

          !2: Analytical solution assuming only \alpha_B and Gamma are constant over update
          else if(ion_type .eq. 2) then
            !Equation 18 of Kim et al. 2017
            xeq_0 = (2*alpha_rec_cell*nH)/(Gamma + 2*alpha_rec_cell*nH + SQRT(Gamma**2 + 4*alpha_rec_cell*nH*Gamma))
            x0    = iony
            xplus = 1./xeq_0
            Kvar = exp(-(xplus - xeq_0)*(dt)*alpha_rec_cell*nH)
            
            if(Gamma .eq. 0) then 
              !If fully neutral, remain neutral
              if(ne .eq. 0) then
                xn = x0 
              else
                !Eq 19 of Kim et al 2017
                xn = (x0+(1-x0)*alpha_rec_cell*nH*(time+dt))/(1+(1-x0)*alpha_rec_cell*nH*(time+dt))
              endif

            else
              ! Eq 17 of Kim et al 2017
              xn = xeq_0 + (xplus - xeq_0)*(x0-xeq_0)*Kvar/((xplus - x0) + (x0- xeq_0)*Kvar)
            endif

          !Implicit Backward-Euler version of Eq. 16 in Kim et al. 2017
          else if(ion_type .eq. 3) then 
            call Implicit_Quadratic_Alt(xn)

          !Implicit Backward-Euler version of above in alternate form that has exact photon conservation
          else if(ion_type .eq. 4) then 
            call Implicit_Quadratic(xn)

          else
            call Driver_abortFlash("[IONIZE] : ionisation type (ion_type) not recognised")

          endif

          !Safety checks (1.e-5 is the tolerance above which I am asssuming something is wrong)
          if(xn .gt. 1 + 1.e-5) then 
           print *, 'xn,xeq_0,x0,dt,tir,tauh=',xn,xeq_0,x0,dt,tir, solnData(TAUH_VAR,i,j,k)
           call Driver_abortFlash("[IONIZE]: Neutral fraction >1")
          endif
          if(xn .lt. -1.e-5) then
            print *, 'xn,xeq_0,x0,dt,tir=',xn,xeq_0,x0,dt,tir
            call Driver_abortFlash("[IONIZE]: Neutral fraction <0")
          endif
          if(isnan(xn)) then
            print *, 'Gamma, IonizationRate, nH0', Gamma, IonizationRate, nH0
            call Driver_abortFlash("[IONIZE]: Nan Neutral fraction")
          endif

          !Now confine the ionisation fraction b/w 0 and 1 (in case it is slightly over)
          xn = min(max(gr_smallx ,xn), 1.0)

          !Save last iteration neutral number density
          solnData(IONP_VAR,i,j,k) = iony

          !Store Solution
#ifdef IHA_SPEC
          !KROME version
          !Update KROME mass fraction of H atoms
          !First record the change in n_Hplus -> this many NEW electrons (per volume) would have formed
          ne_new = (1.-xn) * nH - nHplus
          !Now retreive the new number densities of neutral and ionised hydrogen (keeping in mind nH is same)
          nH0 = xn * nH
          nHplus = (1.-xn) * nH
          ne = ne + ne_new

          !We have to convert neutral fraction to a mass fraction (note nH will remain constant)
          solnData(IHA_SPEC,i,j,k) = nH0 * hA/rho
          solnData(IHP_SPEC,i,j,k) = nHplus * hpA/rho
          solnData(ELEC_SPEC,i,j,k) = ne*elecA/rho

          !Normalise species fraction to maintain closure relation to machine precision
          ! sum up species fractions
          sum = 0.0
          do n = SPECIES_BEGIN, SPECIES_END
            sum = sum + max(gr_smallx, solnData(n,i,j,k))
          enddo

          ! re-normalise sum of species fractions to 1
          do n = SPECIES_BEGIN, SPECIES_END
            solnData(n,i,j,k) = solnData(n,i,j,k) / sum
          enddo
#elif defined(IONY_MSCALAR)
          !OnlyH version: much simpler!
          solnData(IONY_MSCALAR,i,j,k) = xn
#endif

        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

CONTAINS

  !Two implicit versions of the rate equation (NOTE: this is not used by default -- explicit substepping invoked rather.)
  !1. Implicit_Quadratic: Where the ionisation term is not implicitly defined
  !2. Implicit_Quadratic_Alt: The ionisation term is evolving implicitly; more accurate!

  SUBROUTINE Implicit_Quadratic(iony_new)

    real, intent(out) :: iony_new
    real :: a, b, c
    real :: root_neg, root_pos
    real :: ne_curr, ion_curr

    !Current value of ne and neutral fraction
    ne_curr = ne
    ion_curr = iony

    !Coefficients of ax**2 + b*x + c
    !These coefficients are obtained by a Euler backward discretisation
    !ne is set as ne + \delta (nHplus) -- i.e. ne + new electrons -- this makes it fully implicit
    a = alpha_rec_cell * nH * dt
    b = -1*(alpha_rec_cell*ne_curr * dt + (1.+ion_curr)*alpha_rec_cell*nH*dt + 1)
    c = alpha_rec_cell*ne_curr * dt + alpha_rec_cell * nH * ion_curr * dt + iono - IonizationRate * dt/nH

    root_neg = (-b - SQRT(b**2 - 4*a*c))/(2.*a)
    root_pos = (-b + SQRT(b**2 - 4*a*c))/(2.*a)
    
    iony_new = root_neg
    if(iony_new .lt. 0) iony_new = 1.e-6
    != MAX(0.5*iony_new + 0.5 * solnData(IONP_VAR,i,j,k),1.e-6)

    iony_new = iony_new * 0.25 + 0.75*solnData(IONP_VAR,i,j,k)
  END SUBROUTINE Implicit_Quadratic


  SUBROUTINE Implicit_Quadratic_Alt(iony_new)

    real, intent(out) :: iony_new
    real :: a, b, c
    real :: root_neg, root_pos
    real :: ne_curr, ion_curr
    
    !Current value of ne and neutral fraction
    ne_curr = ne
    ion_curr = iony

    !Coefficients of ax**2 + b*x + c
    a = alpha_rec_cell * nH * dt
    b = -1*(alpha_rec_cell*ne_curr * dt + (1.+ion_curr)*alpha_rec_cell*nH*dt + 1 + Gamma*dt)
    c = alpha_rec_cell*ne_curr * dt + alpha_rec_cell * nH * ion_curr * dt + iono

    root_neg = (-b - SQRT(b**2 - 4*a*c))/(2.*a)
    root_pos = (-b + SQRT(b**2 - 4*a*c))/(2.*a)
    
    iony_new = root_neg
    
  END SUBROUTINE Implicit_Quadratic_Alt
END SUBROUTINE IonizeRateEquation

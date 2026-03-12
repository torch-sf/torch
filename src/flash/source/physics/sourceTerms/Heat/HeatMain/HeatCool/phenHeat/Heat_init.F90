!!  ++ / 
!!  + / -   Heating/Cooling
!!   / --
!!written by C. Baczynski, 2012-2013

subroutine Heat_init()

! only few parameters are needed 
! there really should be an except statement
! ah one more: day-umm
  use Heat_data


  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY	: PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, Driver_abortFlash
! use Grid_interface, ONLY:   Grid_setFluxHandling
  use Driver_data, ONLY : dr_dt, dr_restart
  use PhysicalConstants_interface, ONLY :  PhysicalConstants_get

!  use Simulation_data,  ONLY : sim_abar

  use Grid_interface, ONLY : Grid_releaseBlkPtr, Grid_getBlkPtr, Grid_getListOfBlocks

  implicit none

  character(len=80)  :: convert
  integer	:: blockID, thisBlock
	integer :: blockCount
  integer :: blockList(MAXBLOCKS)
  real, pointer, dimension(:,:,:,:)	:: solnData

#include "constants.h"
#include "Flash.h"

! call get_parm_from_context(global_parm_context,"cpnumber", cpnumber)         
! heating and cooling
  call RuntimeParameters_get('tradmin', he_tradmin)
  call RuntimeParameters_get('tradmax', he_tradmax)
  call RuntimeParameters_get('dradmin', he_dradmin)
  call RuntimeParameters_get('dradmax', he_dradmax)
  call RuntimeParameters_get('theatmin', he_theatmin)
  call RuntimeParameters_get('theatmax', he_theatmax)
  call RuntimeParameters_get('absTmin', he_absTmin)
  call RuntimeParameters_get('absTmax', he_absTmax)
  call RuntimeParameters_get('coolOff', he_coolOff)
  call RuntimeParameters_get('h_uv', he_h_uv)
  call RuntimeParameters_get('Gzero', he_Gzero)

! lookit dem constants 
  call PhysicalConstants_get("proton mass", he_protonmass,unitMass="g")
  call PhysicalConstants_get("Boltzmann", he_boltz,unitMass="g",unitLength="cm", & 
					     & unitTime="s")

! switches 
  call RuntimeParameters_get('stratifyHeat', he_stratifyHeat)
  call RuntimeParameters_get('useHeat', he_useHeat)
  call RuntimeParameters_get('subfactor', he_subfactor)
  call RuntimeParameters_get('smallp', he_smallpres)
  call RuntimeParameters_get('dtThres', he_dtThres)

! read in metallicity
  call RuntimeParameters_get('he_abundM', he_abundM)
  call RuntimeParameters_get('he_metal', he_metal)
!  call RuntimeParameters_get('abar', he_abar)

  he_abar = 1.0 + he_abundM*he_metal
  
! CR by JW.  
  call RuntimeParameters_get('use_cr_heating', he_use_cr_heating)
  call RuntimeParameters_get('crIonRate', he_crIonRate)
  call RuntimeParameters_get('crIonNH', he_crIonNH)
  call RuntimeParameters_get('crIonExp', he_crIonExp)
  call RuntimeParameters_get('crIonEnergy', he_crIonEnergy)
  if (he_use_cr_heating) he_crIonEnergy = he_crIonEnergy*1.6022e-12 ! convert energy to cgs.
  
! When to switch off dust heating and cooling stuff. - JW

  call RuntimeParameters_get('dust_sputter_temp', he_dust_sputter_temp)
  call RuntimeParameters_get('he_pe_recipe', he_pe_recipe)
  call RuntimeParameters_get('he_pe_norm', he_pe_norm)
#ifndef IHP_SPEC
  if (he_pe_recipe .eq. "BT94" .or. he_pe_recipe .eq. "WD01") then
    call Driver_abortFlash("[Heat_init] he_pe_recipe=BT94 or WD01 requires ionization from rad xfer unit")
  end if
#endif

  ! Errrybody should know these
  call Driver_getMype(MESH_COMM,he_meshMe)

!  he_abar = sim_abar
  if(he_meshMe==MASTER_PE) print*,"surfin' the cool curve [cooling/heating ON]"

  write(convert, '(i4.4)') he_meshMe
  he_outfile = trim(he_outfile) //trim(convert)// trim('.log')

! init to 0.0

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

#ifdef FERVENT
  do thisBlock = 1, blockCount
    blockID = blockList(thisBlock)
! Get a pointer to solution data 
    call Grid_getBlkPtr(blockID,solnData)
    solnData(PHHE_VAR,:,:,:) = 0d0
    call Grid_releaseBlkPtr(blockID,solnData)
  enddo ! block loop
#endif

   return
end subroutine Heat_init

!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This subroutine performs the radiatiative transfer calculation
!!      for this step using multigroup diffusion theory.
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!   dt     : The time step
!!   pass   : reverses solve direction
!!
!!***

! this corresponds to ionization from FLASH2.5 without the raytracing 
!#define DEBUG_RADTRANS
subroutine RadTrans(nblk, blklst, dt, pass)

! extra unit specific data
  use rt_data, only : rt_protonMass, abu_c, rt_abar, rt_idealgas, rt_dt, rt_dt_pos, rt_maxHchange, &
                      rt_rayTrace

! general radiation transport data
  use RadTrans_data,  ONLY : rt_meshMe
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
      Grid_getBlkIndexLimits, Grid_fillGuardCells

  use Driver_interface, ONLY : Driver_abortFlash
  use Diffuse_interface, ONLY: Diffuse_solveScalar, Diffuse_fluxLimiter

  use RadTrans_data, ONLY: rt_useRadTrans
  use Eos_interface, ONLY: Eos_wrapped
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_rayAdvance
  use Driver_data, ONLY: dr_simTime

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt
  integer, intent(in), optional :: pass

  integer :: j, k, i, l
  real    :: xx, yy, zz
!  solndata
  real, pointer, dimension(:,:,:,:) :: solnData, solnDataCtr
  real, allocatable, dimension(:)   :: xCoord, yCoord, zCoord
  integer                           :: xSizeCoord, ySizeCoord, zSizeCoord
  integer, dimension(2,MDIM)        :: blkLimits, blkLimitsGC
  logical                           :: getGuardCells = .true.

  integer :: tmpID

! for solving 
  real  :: ndens, xH1, xH0, tmp_dt
  real  :: phih, dens, hvphih
  real  :: store, eldens, temp, fac, check
  real  :: hFracNew(0:1), xh(0:1)
  
  !=========================================================================

!  reset radiation timestep
  rt_dt = 1d99

  if(.not. rt_useRadTrans) return

  call Timers_start("RadTrans")

#ifdef DEBUG_RADTRANS
  print*,'entering raytracing'
#endif

! raytracing, should go in Particles_advance.F90 but then
! order in Driver is screwed up, as source term is calculated after call to hydro solver
! as long as the raytracing manipulates the data structure orderly then nothing should happen

  if(rt_rayTrace ) then
    call Timers_start("raytracing")
    call Particles_rayAdvance(dt)
    call Timers_stop("raytracing")
  endif

!  reset radiation timestep
  rt_dt = 1d99

#ifdef DEBUG_RADTRANS
  print*,'leaving raytracing'
#endif

! if single source solving stop here

#ifdef DEBUG_RADTRANS
  print*,'entering rad. solver'
#endif


  call Timers_start("solving") 
!===========================================
! actual radiation transport
!===========================================
  do l = 1, nblk
    tmpID = blklst(l)

! allocate space for dimensions
    call Grid_getBlkPtr(tmpID,solnData)
!		call Grid_getBlkPtr(tmpID,solnDataCtr,SCRATCH_CTR)
    call Grid_getBlkIndexLimits(tmpID,blkLimits,blkLimitsGC)

    xSizeCoord = blkLimitsGC(HIGH,IAXIS)
    ySizeCoord = blkLimitsGC(HIGH,JAXIS)
    zSizeCoord = blkLimitsGC(HIGH,KAXIS)

    allocate(xCoord(xSizeCoord))
    allocate(yCoord(ySizeCoord))
    allocate(zCoord(zSizeCoord))

    call Grid_getCellCoords(IAXIS,tmpID,CENTER,getGuardCells,xCoord,xSizeCoord)
    call Grid_getCellCoords(JAXIS,tmpID,CENTER,getGuardCells,yCoord,ySizeCoord)
    call Grid_getCellCoords(KAXIS,tmpID,CENTER,getGuardCells,zCoord,zSizeCoord)    

  ! loop over all zones in block
    ! for 2d k is just 1
    do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      zz = zCoord(k)
      do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        yy = yCoord(j)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
          xx = xCoord(i)

! not sure if I need the blk pointers or just a solnData array
! temp and density calls for multiple sources here
! additional state variables
          dens = solnData(DENS_VAR,i,j,k)
          temp = solnData(TEMP_VAR,i,j,k)

! neutral and ionised hydrogen fraction
          xH0  = solnData(IHA_SPEC,i,j,k)
          xH1  = solnData(IHP_SPEC,i,j,k)

! phot ionization/ heating rates
          phih   = solnData(PHIO_VAR,i,j,k)

! Convert mass density to number density
          ndens = dens /(rt_abar* rt_protonMass)
          xh = (/xH0,xH1/)

! Initialize new ionization fractions and new temperature to initial ones.
! Needed for convergence and flip-flop check below.
!    			store = xh(0)
          hFracNew = xh
! abu_c for stabiliy (non zero electron density)
          eldens = ndens * hFracNew(1) + abu_c
! Calculate the new and mean ionization state and the new electron
! density.
! hFracNew is changed by calc_ionization, xh is the initial value 
          call calc_ionization(dt, temp, eldens, ndens, hFracNew, xh, phih )

          if(hFracNew(0) + hFracNew(1) .gt. 1.01) then
            print*,phih,dens,temp,hFracNew,xh
            print*,'ion wrong'
            stop ! Don't stop believin... - JW
          endif

! update the ionisation state
          solnData(IHA_SPEC,i,j,k) = hFracNew(0)
          solnData(IHP_SPEC,i,j,k) = hFracNew(1)
          solnData(PHIO_VAR,i,j,k) = 0d0

! reset/treated in phen heating
! 					solnData(PHHE_VAR,i,j,k) = 0d0
			
! time step criterion change of ionized hydrogen fraction in timestep
! this is after changing the ionisation fraction, corrective method
! if change in neutral fraction in current timestep was greater than 0.1 than 
! slow down timestep
          if(rt_rayTrace ) then
            tmp_dt = abs(xH0 - hFracNew(0))
            if(  tmp_dt .gt. rt_maxHchange  ) then
! change in current timstep, next timestep should be smaller 
              fac = rt_maxHchange/tmp_dt
              rt_dt = fac*dt
              rt_dt_pos(1) = 1
              rt_dt_pos(2) = 1
              rt_dt_pos(3) = 1
              rt_dt_pos(4) = tmpID
              rt_dt_pos(5) = rt_meshMe
            endif
          endif
        enddo ! coord loops
      enddo
    enddo

!===========================================
! call to EOS to set zone
!===========================================
! no call to eos hydrogen does not partake in hydrodynamics it is too good for that.
    call Eos_wrapped(MODE_DENS_EI, blkLimits, tmpID)

!  clean up memory 
    call Grid_releaseBlkPtr(tmpID,solnData)
    deallocate(xCoord)
    deallocate(yCoord)
    deallocate(zCoord)
  enddo ! block

#ifdef DEBUG_RADTRANS
	print*,'leaving rad. solver'
#endif

  call Timers_stop("solving") 
  call Timers_stop("RadTrans")

  return
end subroutine RadTrans

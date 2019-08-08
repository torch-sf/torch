!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkAccelGasOnSinks
!!
!! NAME
!!
!!  Particles_sinkAccelGasOnSinks
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelGasOnSinks(OPTIONAL,integer(IN) :: accelProps(MDIM))
!!
!! DESCRIPTION
!!
!!  Computes gas -> sinks gravitational accelerations by direct summation
!!  over all sink particles and grid cells.
!!  For cosmology, will also want to get contribution from PDE
!!  (mapped DM delegate particle density).
!!
!! ARGUMENTS
!!
!!  accelProps : optionally give the indices of the sink particle properties
!!               into which the gas-on-sink accelerations should be stored.
!!               Default is ACCX_PART_PROP, ACCY_PART_PROP, ACCZ_PART_PROP.
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged and renamed to reflect symmetry with Particles_sinkAccelSinksOnGas (Christoph Federrath, 2013)
!!   added optional argument for accel particle properties - Klaus Weide, 2014
!!
!! If accelProps is given but contains an invalid index (e.g., 0), the routine
!! returns without updating any accelerations, but pt_sinkGatherGlobal will
!! still have been called.
!!
!!***

subroutine Gravity_getAccelAtPoint(x, y, z, gaccel_x, gaccel_y, gaccel_z)

!subroutine Particles_sinkAccelGasOnSinks(accelProps)

  use Particles_sinkData
  use pt_sinkSort
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal, pt_sinkEwaldCorrection, &
                              pt_sinkCorrectForPeriodicBCs
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use Driver_data, ONLY : dr_globalMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use tree, ONLY : lnblocks, nodetype
  use Grid_interface, ONLY :  Grid_getCellCoords, Grid_getBlkPhysicalSize, &
                              Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  include "Flash_mpi.h"

  !integer, intent(in), OPTIONAL :: accelProps(MDIM)

  logical, save      :: first_call = .true.
  real, save         :: softening_radius
  real               :: slope, hinv, h2inv, softening_radius_comoving
  real, save         :: newton
  integer            :: i, j, k, p, lb, ierr
  integer            :: size_x, size_y, size_z
  integer, save      :: softeningtype
  real               :: dx_block, dy_block, dz_block, dVol
  real               :: dx, dy, dz, radius, q, kernelvalue, r3, ax, ay, az
  real               :: exc, eyc, ezc
  real               :: prefactor, redshift, oneplusz3
  real               :: size(3), force_sum_x, force_sum_y, force_sum_z, simTime
  character(len=80), save :: softening_type_gas, grav_boundary_type

!  integer, allocatable, dimension(:) :: id_sorted, QSindex
!  real, allocatable, dimension(:) :: ax_sorted, ay_sorted, az_sorted, ax_total, ay_total, az_total
  real,pointer, dimension(:,:,:,: ) :: solnData
  real, dimension(:), allocatable :: xc, yc, zc
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  integer, parameter :: gather_nprops = 5
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP, TAG_PART_PROP, MASS_PART_PROP /)
  integer :: accxProp, accyProp, acczProp

  logical, parameter :: Debug = .false.
  
! Variables for GetGravityAtPoint - Josh Wall 

 real, intent(out) :: gaccel_x, gaccel_y, gaccel_z ! G accelerations to be returned.
 real, intent(in)  :: x, y, z                            ! Location to get gravity at.
 real              :: localgx, localgy, localgz
 
  if (first_call) then

     call RuntimeParameters_get("sink_softening_radius", softening_radius)
     call RuntimeParameters_get("sink_softening_type_gas", softening_type_gas)
     select case (softening_type_gas)
     case ("spline")
        softeningtype=1
     case ("linear")
        softeningtype=2
     case default
        softening_type_gas = "linear"
        softeningtype = 2
        if (dr_globalMe .eq. MASTER_PE) print*, "invalid grav softening type specified"
     end select
     if(dr_globalMe .eq. MASTER_PE) print*, "Particles_sinkKickGas: grav softening type = ", trim(softening_type_gas)

     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

     if ((grav_boundary_type.ne."isolated").and.(grav_boundary_type.ne."periodic")) then
        call Driver_abortFlash("Sink particles can only be used with periodic of isolated gravity type!")
     end if

     call PhysicalConstants_get("Newton", newton)

     first_call = .false.

  end if

  if (Debug .and. dr_globalMe .eq. MASTER_PE) print *, 'Gravity_getAccelAtPoint: entering.'

  call Cosmology_getRedshift(redshift)
  softening_radius_comoving = softening_radius * (1.0 + redshift)

  hinv = 2.0 / softening_radius_comoving
  h2inv = hinv**2
  slope = 1.0 / softening_radius_comoving**3

  oneplusz3 = (1.0 + redshift)**3.0

  ! Exchange particle information
  call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

!  if (present(accelProps)) then
!     accxProp = accelProps(1)
!     accyProp = accelProps(2)
!     acczProp = accelProps(3)
!  else
     accxProp = ACCX_PART_PROP
     accyProp = ACCY_PART_PROP
     acczProp = ACCZ_PART_PROP
!  end if

!  if (accxProp.LE.0 .OR. accyProp.LE.0 .OR. acczProp.LE.0) return

  ! Clear global accelerations
  
  gaccel_x = 0.0
  gaccel_y = 0.0
  gaccel_z = 0.0
  
  localgx = 0.0
  localgy = 0.0
  localgz = 0.0

  ! Loop over blocks
  do lb = 1, lnblocks

     ! only leaf blocks
     if (nodetype(lb) .eq. 1) then

        call Grid_getBlkPtr(lb,solnData)

        call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

        allocate(xc(size_x))
        allocate(yc(size_y))
        allocate(zc(size_z))

        call Grid_getCellCoords(IAXIS, lb, CENTER, .true., xc, size_x)
        call Grid_getCellCoords(JAXIS, lb, CENTER, .true., yc, size_y)
        call Grid_getCellCoords(KAXIS, lb, CENTER, .true., zc, size_z)

        call Grid_getBlkPhysicalSize(lb,size)
        dx_block = size(1)/real(NXB)
        dy_block = size(2)/real(NYB)
        dz_block = size(3)/real(NZB)
        dVol = dx_block*dy_block*dz_block

        ! loop over cells (exclude guard cells)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 prefactor = -newton * solnData(DENS_VAR,i,j,k)  * dVol
!#ifdef PDE_VAR  ! Make sure this isn't called for AMUSE bridging for now. - JW
!                 prefactor = -newton * (solnData(DENS_VAR,i,j,k) + solnData(PDE_VAR,i,j,k)) * dVol
!#endif
                 ! factor of (1+z)^3 needed in cosmological settings:
                 prefactor = prefactor * oneplusz3

                 ! Loop over all particles, local and global
                 !do p = 1, localnpf

                       ! compute relative distances
                       dx = x - xc(i)
                       dy = y - yc(j)
                       dz = z - zc(k)

                       if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)

                       radius = sqrt(dx**2 + dy**2 + dz**2)

                       ! compute accel
                       if (radius .lt. softening_radius_comoving) then
                          if (softeningtype .eq. 1) then    ! spline softening
                             q = radius*hinv
                             if ((q.gt.1.0e-5) .and. (q.lt.1.0)) &
                                & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                             if ((q.ge.1.0)    .and. (q.lt.2.0)) &
                                & kernelvalue = h2inv * &
                                & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                             ax = kernelvalue*dx
                             ay = kernelvalue*dy
                             az = kernelvalue*dz
                          end if

                          if (softeningtype .eq. 2) then ! linear kernel inside softening_radius
                             ax = dx*slope
                             ay = dy*slope
                             az = dz*slope
                          end if
                       else
                          r3 = 1.0 / radius**3
                          ax = dx*r3
                          ay = dy*r3
                          az = dz*r3
                       end if

                       if (grav_boundary_type .eq. "periodic") then
                          call pt_sinkEwaldCorrection(abs(dx), abs(dy), abs(dz), exc, eyc, ezc)
                          ax = ax - sign(exc,dx)
                          ay = ay - sign(eyc,dy)
                          az = az - sign(ezc,dz)
                       endif

                       ! add cell contribution to particle accel
                       localgx = localgx + prefactor*ax
                       localgy = localgy + prefactor*ay
                       localgz = localgz + prefactor*az

                 !end do ! loop over all particles

              enddo  ! i
           enddo  ! j
        enddo  ! k

        call Grid_releaseBlkPtr(lb,solnData)

        deallocate(xc)
        deallocate(yc)
        deallocate(zc)

     end if   ! nodetype

  enddo  ! loop over blocks


  ! allocate temporary arrays

!  allocate(id_sorted(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate id_sorted")
!  allocate(QSindex(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate QSindex")
!  allocate(ax_sorted(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate ax_sorted")
!  allocate(ay_sorted(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate ay_sorted")
!  allocate(az_sorted(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate az_sorted")
!  allocate(ax_total(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate ax_total")
!  allocate(ay_total(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate ay_total")
!  allocate(az_total(localnpf), stat=ierr)
!  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinks:  could not allocate az_total")

!  ! sort global particles list before global all sum
!  do p = 1, localnpf
!     id_sorted(p) = int(particles_global(TAG_PART_PROP,p))
!  enddo

!  call NewQsort_IN(id_sorted, QSindex)

!  ! now particles are sorted by their tag
!  do p = 1, localnpf
!     ax_sorted(p) = particles_global(accxProp, QSindex(p))
!     ay_sorted(p) = particles_global(accyProp, QSindex(p))
!     az_sorted(p) = particles_global(acczProp, QSindex(p))
!     ax_total(p) = 0.0
!     ay_total(p) = 0.0
!     az_total(p) = 0.0
!  enddo

  ! Communicate to get total contribution from all cells on all procs
  call MPI_ALLREDUCE(localgx, gaccel_x, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(localgy, gaccel_y, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(localgz, gaccel_z, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  
!  do p = 1, localnpf
!     particles_global(accxProp, QSindex(p)) = ax_total(p)
!     particles_global(accyProp, QSindex(p)) = ay_total(p)
!     particles_global(acczProp, QSindex(p)) = az_total(p)
!  end do

  ! compute the total force from the gas on all sink particles (only for debugging purposes)
!  if (Debug) then

!      force_sum_x = 0.
!      force_sum_y = 0.
!      force_sum_z = 0.

!      do p = 1, localnpf
!        force_sum_x = force_sum_x + particles_global(MASS_PART_PROP,p)*gaccel_x(p)  !particles_global(accxProp,p)
!        force_sum_y = force_sum_y + particles_global(MASS_PART_PROP,p)*gaccel_y(p)  !particles_global(accyProp,p)
!        force_sum_z = force_sum_z + particles_global(MASS_PART_PROP,p)*gaccel_z(p)  !particles_global(acczProp,p)
!      end do

!      call Driver_getSimTime(simTime)
!      if (dr_globalMe .eq. MASTER_PE) then
!        write(*,'(A,4(1X,E17.10))') &
!          & 'Particles_sinkAccelGasOnSinks: Total force GAS->SINKS (time, x,y,z) = ', &
!          & simTime, force_sum_x, force_sum_y, force_sum_z
!      endif

!  endif

!  do p = 1, localnp
!     particles_local(accxProp,p) = particles_global(accxProp,p)
!     particles_local(accyProp,p) = particles_global(accyProp,p)
!     particles_local(acczProp,p) = particles_global(acczProp,p)
!  end do

!  deallocate(id_sorted)
!  deallocate(QSindex)
!  deallocate(ax_sorted)
!  deallocate(ay_sorted)
!  deallocate(az_sorted)
!  deallocate(ax_total)
!  deallocate(ay_total)
!  deallocate(az_total)

!  if (Debug .and. dr_globalMe .eq. MASTER_PE) print *, 'Particles_sinkAccelGasOnSinks: exiting.'

  return

end subroutine Gravity_getAccelAtPoint

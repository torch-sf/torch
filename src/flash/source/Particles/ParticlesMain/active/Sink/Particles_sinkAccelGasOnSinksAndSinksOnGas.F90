!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkAccelGasOnSinksAndSinksOnGas
!!
!! NAME
!!
!!  Particles_sinkAccelGasOnSinksAndSinksOnGas
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelGasOnSinksAndSinksOnGas(integer, OPTIONAL, intent(IN) :: accelProps(MDIM),
!!                                                  integer, OPTIONAL ,intent(IN) :: accelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  Computes gas -> sinks and sinks -> gas gravitational accelerations
!!  by direct summation over all sink particles and grid cells.
!!  For cosmology, will also want to get contribution from PDE
!!  (mapped DM delegate particle density).
!!
!! ARGUMENTS
!!
!!  accelProps : optionally give the indices of the sink particle properties
!!               into which the gas-on-sink accelerations should be stored.
!!               Default is ACCX_PART_PROP, ACCY_PART_PROP, ACCZ_PART_PROP.
!!
!!   accelVars - optionally give the indices of the UNK variables
!!               into which the sink-on-gas accelerations should be stored.
!!               Default is SGAX_VAR,SGAY_VAR,SGAZ_VAR.
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2015
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged and renamed to reflect symmetry with Particles_sinkAccelSinksOnGas (Christoph Federrath, 2013)
!!   added optional argument for accel particle properties - Klaus Weide, 2014
!!   merged with Particles_sinkAccelSinksOnGas to speed up computation (Christoph Federrath, 2015)
!!
!!
!! If accelProps is given but contains an invalid index (e.g., 0), the routine
!! returns without updating any sink particle accelerations, but gas solution variables
!! will still be updated for the sink-on-gas accelerations (and pt_sinkGatherGlobal will
!! still have been called).
!!
!!***

subroutine Particles_sinkAccelGasOnSinksAndSinksOnGas(accelProps,accelVars)

  use Particles_data,   ONLY: pt_meshComm, pt_globalMe
  use Particles_sinkData,ONLY:particles_global, particles_local, localnp, localnpf, &
                              useSinkParticles
  use pt_sinkSort,      ONLY: NewQsort_IN
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal, pt_sinkEwaldCorrection, &
                              pt_sinkCorrectForPeriodicBCs
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY :  Grid_getCellCoords, Grid_getBlkPhysicalSize, &
                              Grid_getLocalNumBlks, Grid_getBlkType, &
                              Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  include "Flash_mpi.h"

  integer, intent(in), OPTIONAL :: accelProps(MDIM)
  integer, intent(in), OPTIONAL :: accelVars(MDIM)

  logical, save      :: first_call = .true.
  real, save         :: softening_radius
  real               :: slope, hinv, h2inv, softening_radius_comoving
  real               :: pmass, paccx, paccy, paccz
  real, save         :: newton
  integer            :: i, j, k, p, lb, ierr
  integer            :: size_x, size_y, size_z
  integer, save      :: softeningtype
  real               :: dx_block, dy_block, dz_block, dVol
  real               :: dx, dy, dz, radius, q, kernelvalue, r3, ax, ay, az
  real               :: exc, eyc, ezc
  real               :: prefactor_gos, prefactor_sog, redshift, oneplusz3
  real               :: size(3)
  character(len=80), save :: softening_type_gas, grav_boundary_type

  integer*8, allocatable, dimension(:) :: QSindex, id_sorted
  real, allocatable, dimension(:) :: ax_sorted, ay_sorted, az_sorted, ax_total, ay_total, az_total
  real,pointer, dimension(:,:,:,: ) :: solnData
  real, dimension(:), allocatable :: xc, yc, zc
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer                          :: blockCount, blockID
  integer,dimension(:),allocatable :: blockList
  integer :: accxProp, accyProp, acczProp, accxVar, accyVar, acczVar
  logical :: updateSinkProps

  integer, parameter :: gather_nprops = 5
  integer, dimension(gather_nprops), save :: gather_propinds = &
    (/ integer :: POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP, TAG_PART_PROP, MASS_PART_PROP /)

  logical, parameter :: Debug = .false.

  if (first_call) then

    call PhysicalConstants_get("Newton", newton)

    if (useSinkParticles) then
       call RuntimeParameters_get("sink_softening_radius", softening_radius)
       call RuntimeParameters_get("sink_softening_type_gas", softening_type_gas)
       select case (softening_type_gas)
          case ("spline")
             softeningtype = 1
          case ("linear")
             softeningtype = 2
          case default
             softening_type_gas = "linear"
             softeningtype = 2
             if (pt_globalMe .eq. MASTER_PE) print*, 'invalid sink_softening_type_gas specified. using default-> ', &
                                               & trim(softening_type_gas)
       end select
       if (pt_globalMe .eq. MASTER_PE) &
        & print*, 'Particles_sinkAccelSinksOnGasAndSinksOnGas:: sink_softening_type_gas = ', trim(softening_type_gas)

       call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

       if ((grav_boundary_type .ne. "isolated") .and. (grav_boundary_type .ne. "periodic")) then
         call Driver_abortFlash('Sink particles only work with isolated or periodic gravity boundaries.')
       endif

    else
       softening_radius = 0.e0
       slope  = 0.e0
    endif

    first_call = .false.
  endif

  !==============================================================================

  if (localnpf .eq. 0) return

  call Timers_start("AccelGasSinks-SinksGas")

  call Timers_start("loop")

  if (Debug .and. pt_globalMe .eq. MASTER_PE) print *, 'Particles_sinkAccelGasOnSinksAndSinksOnGas: entering.'

  call Cosmology_getRedshift(redshift)
  softening_radius_comoving = softening_radius * (1.0 + redshift)
  hinv  = 2.0/softening_radius_comoving !!! makes sure that only for r < r_soft actual softening occurs
  h2inv = hinv**2
  slope = 1.0 / softening_radius_comoving**3
  oneplusz3 = (1.0 + redshift)**3.0

  prefactor_sog = -newton*oneplusz3

  ! Exchange particle information
  call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

  if (present(accelProps)) then
     accxProp = accelProps(1)
     accyProp = accelProps(2)
     acczProp = accelProps(3)
  else
     accxProp = ACCX_PART_PROP
     accyProp = ACCY_PART_PROP
     acczProp = ACCZ_PART_PROP
  end if
  if (accxProp.le.0 .or. accyProp.le.0 .or. acczProp.le.0) then
     updateSinkProps = .FALSE.
!!$     call Driver_abortFlash('Particles_sinkAccelGasOnSinksAndSinksOnGas: ERROR in index of particle property.')
  else
     updateSinkProps = .TRUE.
  end if

  if (present(accelVars)) then
     accxVar = accelVars(1)
     accyVar = accelVars(2)
     acczVar = accelVars(3)
  else
     accxVar = SGAX_VAR
     accyVar = SGAY_VAR
     acczVar = SGAZ_VAR
  end if
  if (accxVar.le.0 .or. accyVar.le.0 .or. acczVar.le.0) then
     call Driver_abortFlash('Particles_sinkAccelGasOnSinksAndSinksOnGas: ERROR in index of grid variable.')
  end if

  call Grid_notifySolnDataUpdate( (/accxVar,accyVar,acczVar/) )

  ! Clear global accelerations
  if (updateSinkProps) then
     particles_global(accxProp, 1:localnpf) = 0.0
     particles_global(accyProp, 1:localnpf) = 0.0
     particles_global(acczProp, 1:localnpf) = 0.0
  end if

  call Grid_getLocalNumBlks(blockCount)
  allocate(blockList(blockCount))
  call Grid_getListOfBlocks(LEAF, blockList, blockCount)
  ! Loop over leaf blocks
  do lb = 1, blockCount

     blockID = blockList(lb)

        call Grid_getBlkPtr(blockID,solnData)

        call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

        allocate(xc(size_x))
        allocate(yc(size_y))
        allocate(zc(size_z))

        call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xc, size_x)
        call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yc, size_y)
        call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zc, size_z)

        call Grid_getBlkPhysicalSize(blockID,size)
        dx_block = size(1)/real(NXB)
        dy_block = size(2)/real(NYB)
        dz_block = size(3)/real(NZB)
        dVol = dx_block*dy_block*dz_block

        ! loop over cells (exclude guard cells)
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 prefactor_gos = -newton * solnData(DENS_VAR,i,j,k)  * dVol
#ifdef PDE_VAR
                 prefactor_gos = -newton * (solnData(DENS_VAR,i,j,k) + solnData(PDE_VAR,i,j,k)) * dVol
#endif
                 ! factor of (1+z)^3 needed in cosmological settings:
                 prefactor_gos = prefactor_gos * oneplusz3

                 ! reset particle accelerations
                 paccx = 0.0
                 paccy = 0.0
                 paccz = 0.0

                 ! Loop over all particles, local and global
                 do p = 1, localnpf

                       ! particle mass
                       pmass = particles_global(MASS_PART_PROP,p)

                       ! compute relative distances
                       dx = particles_global(POSX_PART_PROP,p) - xc(i)
                       dy = particles_global(POSY_PART_PROP,p) - yc(j)
                       dz = particles_global(POSZ_PART_PROP,p) - zc(k)

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

                       if (updateSinkProps) then
                       ! gas on sinks: add cell contribution to particle acceleration
                          particles_global(accxProp,p) = particles_global(accxProp,p) + & 
                            prefactor_gos * ax
                          particles_global(accyProp,p) = particles_global(accyProp,p) + & 
                            prefactor_gos * ay
                          particles_global(acczProp,p) = particles_global(acczProp,p) + & 
                            prefactor_gos * az
                       end if

                       ! sinks on gas accelerations
                       paccx = paccx - ax * pmass
                       paccy = paccy - ay * pmass
                       paccz = paccz - az * pmass

                 end do ! loop over all particles

                 ! sinks on gas: x-acceleration:
                 solnData(accxVar,i,j,k) = paccx * prefactor_sog
                 ! sinks on gas: y-acceleration:
                 solnData(accyVar,i,j,k) = paccy * prefactor_sog
                 ! sinks on gas: z-acceleration:
                 solnData(acczVar,i,j,k) = paccz * prefactor_sog

              enddo  ! i
           enddo  ! j
        enddo  ! k

        call Grid_releaseBlkPtr(blockID,solnData)

        deallocate(xc)
        deallocate(yc)
        deallocate(zc)


  enddo  ! loop over blocks

  deallocate(blockList)

  call Timers_stop("loop")

  if (updateSinkProps) then

  call Timers_start("reduce")

  ! allocate temporary arrays

  allocate(id_sorted(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate id_sorted")
  allocate(QSindex(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate QSindex")
  allocate(ax_sorted(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate ax_sorted")
  allocate(ay_sorted(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate ay_sorted")
  allocate(az_sorted(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate az_sorted")
  allocate(ax_total(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate ax_total")
  allocate(ay_total(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate ay_total")
  allocate(az_total(localnpf), stat=ierr)
  if (ierr.ne.0) call Driver_abortFlash ("Particles_sinkAccelGasOnSinksAndSinksOnGas:  could not allocate az_total")

  ! sort global particles list before global all sum
  do p = 1, localnpf
     id_sorted(p) = int(particles_global(TAG_PART_PROP,p),8)
  enddo

  call NewQsort_IN(id_sorted, QSindex)

  ! now particles are sorted by their tag
  do p = 1, localnpf
     ax_sorted(p) = particles_global(accxProp, QSindex(p))
     ay_sorted(p) = particles_global(accyProp, QSindex(p))
     az_sorted(p) = particles_global(acczProp, QSindex(p))
     ax_total(p) = 0.0
     ay_total(p) = 0.0
     az_total(p) = 0.0
  enddo

  ! Communicate to get total contribution from all cells on all procs
  call MPI_ALLREDUCE(ax_sorted, ax_total, localnpf, FLASH_REAL, MPI_SUM, pt_meshComm, ierr)
  call MPI_ALLREDUCE(ay_sorted, ay_total, localnpf, FLASH_REAL, MPI_SUM, pt_meshComm, ierr)
  call MPI_ALLREDUCE(az_sorted, az_total, localnpf, FLASH_REAL, MPI_SUM, pt_meshComm, ierr)

  do p = 1, localnpf
     particles_global(accxProp, QSindex(p)) = ax_total(p)
     particles_global(accyProp, QSindex(p)) = ay_total(p)
     particles_global(acczProp, QSindex(p)) = az_total(p)
  end do

  do p = 1, localnp
     particles_local(accxProp,p) = particles_global(accxProp,p)
     particles_local(accyProp,p) = particles_global(accyProp,p)
     particles_local(acczProp,p) = particles_global(acczProp,p)
  end do

  deallocate(id_sorted)
  deallocate(QSindex)
  deallocate(ax_sorted)
  deallocate(ay_sorted)
  deallocate(az_sorted)
  deallocate(ax_total)
  deallocate(ay_total)
  deallocate(az_total)

  if (Debug .and. pt_globalMe .eq. MASTER_PE) print *, 'Particles_sinkAccelGasOnSinksAndSinksOnGas: exiting.'

  call Timers_stop("reduce")

  end if

  call Timers_stop("AccelGasSinks-SinksGas")

  return

end subroutine Particles_sinkAccelGasOnSinksAndSinksOnGas

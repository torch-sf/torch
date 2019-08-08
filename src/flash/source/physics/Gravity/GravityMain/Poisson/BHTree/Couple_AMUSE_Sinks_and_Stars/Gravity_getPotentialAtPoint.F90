!!****if* source/Gravity/GravityMain/Poisson/Multigrid
!!
!! NAME
!!
!!  Gravity_getPotentialAtPoint (Multigrid)
!!
!! SYNOPSIS
!!
!!  call Gravity_getPotentialAtPoint
!!
!! DESCRIPTION
!!
!!  Computes the gravitational potential from the gas using
!!  the tree implementation by K. Olson. Based on the file
!!  Particles_getAccelGasOnSinks by R. Banerjee.
!!
!! ARGUMENTS
!!
!! x, y, z : Location at which to calculate the gravity.
!! grav_pot : Returned gravitational potential at x,y,z.
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged and renamed to reflect symmetry with Particles_sinkAccelSinksOnGas (Christoph Federrath, 2013)
!!   adapted for current use by J. Wall 2015
!!***

subroutine Gravity_getPotentialAtPoint(x, y, z, grav_pot_total)

  use Particles_sinkData
  use pt_sinkSort
  use pt_sinkInterface, only: pt_sinkGatherGlobal
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

  logical, save      :: first_call = .true.
  real, save         :: softening_radius
  real               :: slope, hinv, h2inv, softening_radius_comoving
  real, save         :: maxradius_pbc, newton
  integer            :: i, j, k, p, lb, ierr, n, nx, ny, nz
  real, save         :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz
  integer, parameter :: nrep_pbc=6
  real, save         :: LxPBC(-nrep_pbc:nrep_pbc), LyPBC(-nrep_pbc:nrep_pbc), LzPBC(-nrep_pbc:nrep_pbc)
  integer            :: size_x, size_y, size_z
  integer, save      :: softeningtype
  real               :: dx_block, dy_block, dz_block, dVol
  real               :: dx, dy, dz, x2, y2, z2, radius, q, kernelvalue, r3, ax, ay, az
  real               :: dx_inside, dy_inside, dz_inside
  real               :: prefactor, redshift, oneplusz3
  real               :: size(3), force_sum_x, force_sum_y, force_sum_z, simTime
  character(len=80), save :: softening_type_gas, grav_boundary_type

  integer, allocatable, dimension(:) :: id_sorted, QSindex
  real, allocatable, dimension(:) :: ax_sorted, ay_sorted, az_sorted, ax_total, ay_total, az_total
  real,pointer, dimension(:,:,:,: ) :: solnData
  real, dimension(:), allocatable :: xc, yc, zc
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  logical, parameter :: Debug = .false.
  
! For Bridge implementation with AMUSE - Josh Wall
 real    :: pot_kernel, pgpot
 
! Variables for GetGravityAtPoint - Josh Wall

 real     :: grav_pot 
 real     :: grav_pot_total  ! G  pot to be returned.
 real     :: x, y, z         ! Location to get gravity at.

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
     if(dr_globalMe .eq. MASTER_PE) print*, "Particles_sinkAccelGasOnSinks: grav softening type = ", trim(softening_type_gas)

     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

     if ((grav_boundary_type.ne."isolated").and.(grav_boundary_type.ne."periodic")) then
        call Driver_abortFlash("Sink particles can only be used with perioidic of isolated gravity type!")
     end if

     call RuntimeParameters_get("xmin", xmin)
     call RuntimeParameters_get("xmax", xmax)
     call RuntimeParameters_get("ymin", ymin)
     call RuntimeParameters_get("ymax", ymax)
     call RuntimeParameters_get("zmin", zmin)
     call RuntimeParameters_get("zmax", zmax)

     Lx = xmax-xmin
     Ly = ymax-ymin
     Lz = zmax-zmin

     if (grav_boundary_type .eq. "periodic") then
        maxradius_pbc = real(nrep_pbc)*min(Lx,Ly,Lz)
        if(dr_globalMe .eq. MASTER_PE) print*, "Particles_sinkAccelGasOnSinks: maxradius_pbc = ", maxradius_pbc
        do n = -nrep_pbc, nrep_pbc
           LxPBC(n) = n*Lx
           LyPBC(n) = n*Ly
           LzPBC(n) = n*Lz
        enddo
     end if

     call PhysicalConstants_get("Newton", newton)

     first_call = .false.

  end if

  call Cosmology_getRedshift(redshift)
  softening_radius_comoving = softening_radius * (1.0 + redshift)

  hinv = 2.0 / softening_radius_comoving
  h2inv = hinv**2
  slope = 1.0 / softening_radius_comoving**3

  oneplusz3 = (1.0 + redshift)**3.0


     grav_pot_total = 0.0
     grav_pot = 0.0 ! And potentials. -JW


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
#ifdef PDE_VAR
                 prefactor = -newton * (solnData(DENS_VAR,i,j,k) + solnData(PDE_VAR,i,j,k)) * dVol
#endif
                 ! factor of (1+z)^3 needed in cosmological settings:
                 prefactor = prefactor * oneplusz3

                 if (grav_boundary_type .eq. "isolated") then

                    ! Loop over all particles, local and global
                    !do p = 1, localnpf

                       ! compute relative distances
                       dx = x - xc(i)
                       dy = y - yc(j)
                       dz = z - zc(k)
                       radius = sqrt(dx*dx+dy*dy+dz*dz)

                       ! compute pot
                       if (radius .lt. softening_radius_comoving) then
                          if(softeningtype .eq. 1) then    ! spline softening
                             q = radius*hinv
                             if ((q.gt.1.0e-5) .and. (q.lt.1.0)) &
                                & pot_kernel  = hinv*(2.0/3.0*q**2.0 - 1.0/10.0*q**4.0 + &
                                & 1.0/10.0*q**5.0 - 7.0/5.0)
                             if ((q.ge.1.0)    .and. (q.lt.2.0)) &
                                & pot_kernel = hinv*(4.0/3.0*q**2.0 - q**3.0 + 3.0/10.0*q**4.0 - &
                                & 1/30.0*q**5.0 - 8.0/5.0 +1.0/(15.0*q))

                             pgpot = pot_kernel ! Again -JW
                          end if

                          if (softeningtype .eq. 2) then ! linear kernel inside softening_radius
                             pgpot = slope
                          end if
                          
                       else
                          pgpot = 1/radius ! -JW

                       end if

                       ! add cell contribution to pot
                       grav_pot = grav_pot + prefactor*pgpot

                 end if

                 if (grav_boundary_type .eq. "periodic") then


                       pgpot = 0.0

                       ! compute relative distances
                       dx_inside = x - xc(i)
                       dy_inside = y - yc(j)
                       dz_inside = z - zc(k)

                       do nx = -nrep_pbc, nrep_pbc
                          dx = dx_inside + LxPBC(nx)
                          x2 = dx*dx
                          do ny = -nrep_pbc, nrep_pbc
                             dy = dy_inside + LyPBC(ny)
                             y2 = dy*dy
                             do nz = -nrep_pbc, nrep_pbc
                                dz = dz_inside + LzPBC(nz)
                                z2 = dz*dz

                                radius = sqrt(x2 + y2 + z2)

                                if (radius .lt. maxradius_pbc) then

                                   if (radius .lt. softening_radius_comoving) then

                                      if(softeningtype .eq. 1) then    ! spline softening
                                         q = radius*hinv
                                         if ((q.gt.1.0e-5) .and. (q.lt.1.0)) &
                                            & pot_kernel  = hinv*(2.0/3.0*q**2.0 - 1.0/10.0*q**4.0 + &
                                            & 1.0/10.0*q**5.0 - 7.0/5.0)
                                         if ((q.ge.1.0)    .and. (q.lt.2.0)) &
                                            & pot_kernel = hinv*(4.0/3.0*q**2.0 - q**3.0 + 3.0/10.0*q**4.0 - &
                                            & 1/30.0*q**5.0 - 8.0/5.0 +1.0/(15.0*q))

                                         pgpot = pot_kernel ! Again -JW
                                      end if

                                      if (softeningtype .eq. 2) then ! linear kernel inside softening_radius
                                         pgpot = slope
                                      end if
                                   else
                                      pgpot = 1/radius
                                   endif

                                end if ! within maxradius_pbc

                             enddo  ! nz
                          enddo    ! ny
                       enddo      ! nz

                       ! add cell contribution to pot
                       grav_pot = grav_pot + prefactor*pgpot


                 endif    ! gravity boundary

              enddo  ! i
           enddo  ! j
        enddo  ! k

        call Grid_releaseBlkPtr(lb,solnData)

        deallocate(xc)
        deallocate(yc)
        deallocate(zc)

     end if   ! nodetype

  enddo  ! loop over blocks

  ! Communicate to get total contribution from all cells on all procs
  call MPI_ALLREDUCE(grav_pot, grav_pot_total, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  return

end subroutine Gravity_getPotentialAtPoint

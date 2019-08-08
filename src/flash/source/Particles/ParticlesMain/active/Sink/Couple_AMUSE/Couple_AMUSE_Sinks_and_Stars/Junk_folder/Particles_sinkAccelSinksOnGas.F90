!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkAccelSinksOnGas
!!
!! NAME
!!
!!  Particles_sinkAccelSinksOnGas
!!
!! SYNOPSIS
!!
!!  call Particles_sinkAccelSinksOnGas(integer,intent(IN) :: blockcount,
!!              integer,dimension(blockCount),intent(IN)  :: blocklist,
!!              integer, OPTIONAL ,intent(IN)             :: accelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  Computes SGAX,SGAY,SGAZ UNK vars from sink particles. (sinks -> gas accelerations).
!!  if accelVars is given, store into those UNK vars instead.
!!  Computes sinks -> gas gravitational accelerations by direct summation
!!  over all sink particles and grid cells.
!!
!! ARGUMENTS
!!
!!   blockcount - the number of blocks on this processor
!!
!!   blocklist - the list of blocks held by this processor
!!
!!   accelVars - optionally give the indices of the UNK variables
!!               into which the sink-on-gas accelerations should be stored.
!!               Default is SGAX_VAR,SGAY_VAR,SGAZ_VAR.
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   debugged by Christoph Federrath, 2013
!!   added optional argument for accel vars - Klaus Weide, 2014
!!
!!***

subroutine Particles_sinkAccelSinksOnGas(blockCount,blockList, accelVars)

!==============================================================================

 use Particles_sinkData
 use pt_sinkInterface, ONLY: pt_sinkEwaldCorrection, pt_sinkCorrectForPeriodicBCs
 use Driver_data, ONLY : dr_globalMe
 use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
 use Grid_interface, ONLY : Grid_getCellCoords, Grid_getBlkIndexLimits, & 
     Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkPhysicalSize, &
     Grid_notifySolnDataUpdate
 use RuntimeParameters_interface, ONLY : RuntimeParameters_get
 use Cosmology_interface, ONLY : Cosmology_getRedshift
 use PhysicalConstants_interface, ONLY : PhysicalConstants_get

 implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
include "Flash_mpi.h"

 integer,intent(IN) :: blockCount
 integer,dimension(blockCount),intent(IN) :: blockList
 integer, intent(in), OPTIONAL :: accelVars(MDIM)

 real, POINTER, DIMENSION(:,:,:,:) :: solnData

 character(len=80), save :: softening_type_gas, grav_boundary_type

 logical, save :: first_call = .true.
 real, save    :: newton
 integer, save :: softeningtype
 real, save    :: softening_radius
 real          :: slope, hinv, h2inv
 real          :: radius, prefactor, q, kernelvalue, r3
 real          :: paccx, paccy, paccz, mass
 integer       :: nxbBlock, nybBlock, nzbBlock, ierr, lb
 real          :: blockSize(MDIM), dx_block, dy_block, dz_block, dVol
 real          :: dx, dy, dz, ax, ay, az, exc, eyc, ezc
 real          :: force_sum_x, force_sum_y, force_sum_z, density, simTime
 real          :: force_sum_x_glob, force_sum_y_glob, force_sum_z_glob
 integer       :: ii, jj, kk, p
 integer       :: size_x, size_y, size_z
 integer       :: accxVar, accyVar, acczVar
 real          :: redshift, softening_radius_comoving, oneplusz3

 real, dimension(:), allocatable :: x, y, z
 integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

 logical, parameter :: Debug = .false.

 !==============================================================================

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
            if (dr_globalMe .eq. MASTER_PE) print*, 'invalid sink_softening_type_gas specified. using default-> ', &
                                              & trim(softening_type_gas)
      end select
      if (dr_globalMe .eq. MASTER_PE) print*, 'Particles_sinkAccelSinksOnGas:: sink_softening_type_gas = ', &
                                        & trim(softening_type_gas)

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


! This shouldn't be called when we use the AMUSE coupling. -JW
!  print*, "[Particles_sinkAccelSinksOnGas]: Who be calling me?"
!  return

 !==============================================================================

  if (present(accelVars)) then
     accxVar = accelVars(1)
     accyVar = accelVars(2)
     acczVar = accelVars(3)
  else
     accxVar = SGAX_VAR
     accyVar = SGAY_VAR
     acczVar = SGAZ_VAR
  end if

 if (localnpf .EQ. 0) return

 if (Debug .and. dr_globalMe .eq. MASTER_PE) print *, 'Particles_sinkAccelSinksOnGas: entering.'

 call Cosmology_getRedshift(redshift)
 softening_radius_comoving = softening_radius * (1.0 + redshift)
 hinv  = 2.0/softening_radius_comoving !!! makes sure that only for r < r_soft actual softening occurs
 h2inv = hinv**2
 slope = 1.0/softening_radius_comoving**3
 oneplusz3 = (1.0 + redshift)**3.0
 prefactor = -newton*oneplusz3

 if (Debug) then
   force_sum_x = 0.
   force_sum_y = 0.
   force_sum_z = 0.
 endif

 call Grid_notifySolnDataUpdate( (/accxVar,accyVar,acczVar/) )

 ! Loop through blocks
 do lb = 1, blockCount

    call Grid_getBlkPtr(blockList(lb),solnData)

    call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
    call Grid_getBlkPhysicalSize(blockList(lb), blockSize)

    size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
    size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
    size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

    nxbBlock = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
    nybBlock = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
    nzbBlock = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1

    allocate(x(size_x))
    allocate(y(size_y))
    allocate(z(size_z))

    call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., x, size_x)
    call Grid_getCellCoords(JAXIS, blockList(lb), CENTER, .true., y, size_y)
    call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, .true., z, size_z)

    if (Debug) then
      dx_block = blockSize(1)/real(NXB)
      dy_block = blockSize(2)/real(NYB)
      dz_block = blockSize(3)/real(NZB)
      dVol = dx_block*dy_block*dz_block
    endif

    ! Loop through cells, neglecting one layer of guard cells
    do ii = blkLimitsGC(LOW,IAXIS)+1, blkLimitsGC(HIGH,IAXIS)-1
       do jj = blkLimitsGC(LOW,JAXIS)+1, blkLimitsGC(HIGH,JAXIS)-1
          do kk = blkLimitsGC(LOW,KAXIS)+1, blkLimitsGC(HIGH,KAXIS)-1

               paccx = 0.0
               paccy = 0.0
               paccz = 0.0

               ! Loop over global sink particles
               do p = 1, localnpf

                  mass = particles_global(MASS_PART_PROP,p)

                  dx = x(ii) - particles_global(POSX_PART_PROP,p)
                  dy = y(jj) - particles_global(POSY_PART_PROP,p)
                  dz = z(kk) - particles_global(POSZ_PART_PROP,p)

                  if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)

                  radius = sqrt(dx**2 + dy**2 + dz**2)

                  if (radius .lt. softening_radius_comoving) then
                     if (softeningtype.eq.1) then ! spline softening (see e.g., Price & Monaghan 2007)
                        q = radius*hinv
                        if ((q.gt.1.e-5).and.(q.lt.1.0)) &
                           & kernelvalue = h2inv*(4.0/3.0*q-1.2*q**3+0.5*q**4)/radius
                        if ((q.ge.1.0)  .and.(q.lt.2.0)) & 
                           & kernelvalue = h2inv * &
                           & (8.0/3.0*q-3.0*q**2+1.2*q**3-1.0/6.0*q**4-1.0/(15.0*q**2))/radius
                        ax = kernelvalue*dx
                        ay = kernelvalue*dy
                        az = kernelvalue*dz
                     endif
                     if (softeningtype.eq.2) then ! linear kernel inside smoothing radius
                        ax = dx*slope
                        ay = dy*slope
                        az = dz*slope
                     endif
                  else
                     r3 = 1.0/radius**3
                     ax = dx*r3
                     ay = dy*r3
                     az = dz*r3
                  endif

                  if (grav_boundary_type .eq. "periodic") then
                     call pt_sinkEwaldCorrection(abs(dx), abs(dy), abs(dz), exc, eyc, ezc)
                     paccx = paccx + (ax - sign(exc,dx))*mass
                     paccy = paccy + (ay - sign(eyc,dy))*mass
                     paccz = paccz + (az - sign(ezc,dz))*mass
                  else
                     paccx = paccx + ax*mass
                     paccy = paccy + ay*mass
                     paccz = paccz + az*mass
                  endif

               end do   ! loop over all sinks

               ! x-acceleration:
               solnData(accxVar,ii,jj,kk) = paccx*prefactor

               ! y-acceleration:
               solnData(accyVar,ii,jj,kk) = paccy*prefactor

               ! z-acceleration:
               solnData(acczVar,ii,jj,kk) = paccz*prefactor

               ! compute the total force from the sinks on the gas (only for debugging purposes)
               if (Debug) then

                  density = solnData(DENS_VAR,ii,jj,kk)
#ifdef PDE_VAR
                  density = solnData(DENS_VAR,ii,jj,kk) + solnData(PDE_VAR,ii,jj,kk)
#endif
                  if ( (ii.ge.blkLimits(LOW,IAXIS)) .and. (ii.le.blkLimits(HIGH,IAXIS)) .and. &
                       (jj.ge.blkLimits(LOW,JAXIS)) .and. (jj.le.blkLimits(HIGH,JAXIS)) .and. &
                       (kk.ge.blkLimits(LOW,KAXIS)) .and. (kk.le.blkLimits(HIGH,KAXIS)) ) then
                         force_sum_x = force_sum_x + solnData(accxVar,ii,jj,kk)*density*dVol
                         force_sum_y = force_sum_y + solnData(accyVar,ii,jj,kk)*density*dVol
                         force_sum_z = force_sum_z + solnData(acczVar,ii,jj,kk)*density*dVol
                  endif ! in active region

               endif ! Debug

          enddo   ! cells
       enddo   ! cells
    enddo   ! cells

    deallocate(x)
    deallocate(y)
    deallocate(z)

    call Grid_releaseBlkPtr(blockList(lb),solnData)

 end do   ! blocks


 if (Debug) then

     ! Communicate to get total contribution from all cells on all procs
     call MPI_ALLREDUCE(force_sum_x, force_sum_x_glob, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(force_sum_y, force_sum_y_glob, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(force_sum_z, force_sum_z_glob, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call Driver_getSimTime(simTime)
     if (dr_globalMe .eq. MASTER_PE) &
       & write(*,'(A,4(1X,E17.10))') 'Particles_sinkAccelSinksOnGas: Total force SINKS->GAS (time, x,y,z) = ', &
       & simTime, force_sum_x_glob, force_sum_y_glob, force_sum_z_glob

 endif

 if (Debug .and. dr_globalMe .eq. MASTER_PE) print *, 'Particles_sinkAccelSinksOnGas: exiting.'

 return

end subroutine Particles_sinkAccelSinksOnGas

!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkMarkRefineDerefine
!!
!! NAME
!!
!!  Particles_MarkRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Particles_MarkRefineDerefine()
!!
!! DESCRIPTION
!!
!!  This routine takes care of grid refinement for particles which have
!!  feedback (active/massive particles with either FERVENT or my wind
!!  module included). It enforces refinement in the blocks that contain
!!  the particles, and in the case of the winds, in any block where the
!!  wind will be injected. Note here that FERVENT refinement simply
!!  looks to see if the particle has rays, while the winds check if the
!!  particle exceeds the minimum mass for a wind. - JW
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   cleaned by Christoph Federrath, 2013
!!   adapted for active particles J Wall 2016-2017
!!
!!***
!#define debug
subroutine Particles_MarkRefineDerefine()
#include "Flash.h"
#include "Particles.h"
#include "constants.h"
#include "Multispecies.h"

    use RuntimeParameters_interface, only: RuntimeParameters_get
    use Particles_data, only : particles, pt_numLocal, pt_typeInfo
#ifdef WIND_INJ
    use Particles_windData, only : ref_radius, min_wind_mass, min_jet_mass, max_jet_mass 
    !Added jet params -SA 20230918
#endif

    use Particles_interface, only : Particles_getGlobalNum
  
    use tree
    use paramesh_dimensions
    use physicaldata, ONLY : unk
    use Grid_data, ONLY : gr_maxRefine

    use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPhysicalSize, & 
         Grid_getCellCoords, Grid_getBlkIndexLimits, Grid_getBlkCenterCoords
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Driver_data, ONLY : dr_globalMe, dr_globalComm, dr_globalNumProcs
    use Grid_interface, ONLY : Grid_getBlkPhysicalSize

    use pt_sinkInterface, only: pt_sinkCorrectForPeriodicBCs
    
    implicit none

    ! Arguments

  ! Local data

    integer :: b, ii, jj, kk, p, p_blknum

    real, dimension(:), allocatable :: xc, yc, zc
    integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
    integer :: size_x, size_y, size_z, kp, jp, ip
    logical :: p_found

    real, dimension(MDIM) :: blockSize
    real          :: rad, distx, disty, distz, delta

! Local variable for the minimum wind/jet mass -SA 20230918
    real*8   :: loc_feedback_mass

! Global vars for particle locations.
    real*8, allocatable :: x(:), y(:), z(:), &
                       dmdt(:), v_wind(:), c_time(:), bgdy(:)
    integer    :: w_num
! Local vars for particle locatiosn.
    real, allocatable      :: locx(:), locy(:), locz(:)

! For counting particles.
    integer                :: p_begin, p_end, p_num, p_globalnum, w_numloc

! For MPI comm
    integer                :: num_array(dr_globalNumProcs), &
                              disp(dr_globalNumProcs), ierr, i, &
                              rank_minus_one

logical, save              :: first_call = .true.
character(len=80), save    :: grav_boundary_type

#include "Flash_mpi.h"


! If neither feedback method is in, just return.
! If in jets branch WIND_INJ includes both winds and jets -SA 202309118
#if !defined(FERVENT) && !defined(WIND_INJ)
return
#endif

#ifdef WIND_INJ
if (first_call) then
  call RuntimeParameters_get("ref_radius", ref_radius)
  call Grid_getMinCellSize(delta)
  if (ref_radius < 0.0) &
      ref_radius = 10.0*delta !3.5d0*sqrt(3.0d0)*delta
      !Update to larger value needed for jets -SA 20230918
  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
  call RuntimeParameters_get("min_wind_mass", min_wind_mass)
  call RuntimeParameters_get("min_jet_mass", min_jet_mass) ! -SA 20230918
  call RuntimeParameters_get("max_jet_mass", max_jet_mass)
  first_call = .false.
end if

! Also need to account for jet mass: -SA 20230918
! DEfine a local min feedback mass
loc_feedback_mass = MIN(min_wind_mass, min_jet_mass)
#endif

! Local number of massive/active particles.
  p_begin = pt_typeInfo(PART_TYPE_BEGIN,ACTIVE_PART_TYPE)
  p_num   = pt_typeInfo(PART_LOCAL,ACTIVE_PART_TYPE)
  p_end   = p_num + p_begin - 1
#ifdef debug
print*, "p_num =", p_num, dr_globalMe
#endif

call Particles_getGlobalNum(p_globalnum)

allocate(locx(p_globalnum), locy(p_globalnum), locz(p_globalnum))

w_numloc  = 0
w_num     = 0
num_array = 0

locx = 0.0d0; locy=0.0d0; locz=0.0d0

do p = p_begin, p_end
  
  if (particles(MASS_PART_PROP, p) .gt. loc_feedback_mass) then   
     
     w_numloc = w_numloc + 1
     locx(w_numloc)      = particles(POSX_PART_PROP, p)
     locy(w_numloc)      = particles(POSY_PART_PROP, p)
     locz(w_numloc)      = particles(POSZ_PART_PROP, p)
  
  end if

end do

! Now use MPI to vector gather all the information for how to inject
! the winds on each processor.
#ifdef debug
print*, "w_numloc =", w_numloc, dr_globalMe
#endif

disp = 0
rank_minus_one = dr_globalNumProcs - 1

! Gather the array on the root process. Note that we require the
! user to pass the proper length of the final array. This can be
! gotten from get_number_of_new_tags.

! Make an array of the # of incoming particles from each processor.
call MPI_AllGather(w_numloc, 1, MPI_INTEGER, &
	      num_array, 1, MPI_INTEGER, &
	      dr_globalComm, ierr)

! Allocate the actual arrays to pass.

w_num = sum(num_array)

#ifdef debug
print*, "w_num =", w_num, dr_globalMe
#endif

if (allocated(x)) &
    deallocate(x, y, z)

allocate(x(w_num), y(w_num), z(w_num))

x=0.0d0; y=0.0d0; z=0.0d0

! Set the displacement for the incoming data based on how many
! particles are coming in from each processor. Note the displacement
! for the root process is zero, for rank 1 disp = num on root,
! for rank 2 disp = num on root + num on 1, etc etc.

do i=1, dr_globalNumProcs-1

  disp(i+1) = disp(i) + num_array(i)

end do
#ifdef debug
print*, "About to gather.", dr_globalMe
print*, "num_array =", num_array, dr_globalMe
print*, "disp =", disp, dr_globalMe
#endif
! Now actually gather the info on each proc using the variable length array
! gather command in MPI.
call MPI_AllGatherv(locx, w_numloc, FLASH_REAL, x, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locy, w_numloc, FLASH_REAL, y, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
call MPI_AllGatherv(locz, w_numloc, FLASH_REAL, z, num_array, &
	       disp, FLASH_REAL, dr_globalComm, ierr)
         
#ifdef debug
print*, "x =", x
print*, "y =", y
print*, "z =", z
#endif

#ifdef FERVENT
       ! Any block with a feedback particle in it should be at highest refinement level. - JW
       ! This part is basically for HII regions to be well resolved initially
       ! for radiation feedback.
       do p = 1, pt_numLocal
          if (particles(EION_PART_PROP, p) > 0.0d0) then 
            p_blknum = int(particles(BLK_PART_PROP, p))
            if (lrefine(p_blknum) .lt. gr_maxRefine) then
               refine(p_blknum) = .true.
               derefine(p_blknum) = .false.
               stay(p_blknum) = .true.

            end if
            if (lrefine(p_blknum) .eq. gr_maxRefine) then
               derefine(p_blknum) = .false.
               stay(p_blknum) = .true.
            end if
          end if
       end do
#endif
       ! Any cell within accretion_radius of sink particle should be at the
       ! highest refinement level (its block, to be precise)
#if defined(WIND_INJ)
       do b = 1, lnblocks
          ! cell within injection radius?
          p_found = .false.
          
          if (nodetype(b).eq.1) then

             ! find cells (including GCs) within sink particle accretion radius

             call Grid_getBlkIndexLimits(b, blkLimits, blkLimitsGC)
             size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
             size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
             size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

             allocate(xc(size_x))
             allocate(yc(size_y))
             allocate(zc(size_z))

             call Grid_getCellCoords(IAXIS, b, CENTER, .true., xc, size_x)
             call Grid_getCellCoords(JAXIS, b, CENTER, .true., yc, size_y)
             call Grid_getCellCoords(KAXIS, b, CENTER, .true., zc, size_z)

             do kp = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
                do jp = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
                   do ip = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

                      ! For massive particles globally across simulation.  
                      do p = 1, w_num
                         distx = xc(ip) - x(p)
                         disty = yc(jp) - y(p)
                         distz = zc(kp) - z(p)
                         if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
                         rad = sqrt(distx**2 + disty**2 + distz**2)

                         ! Added for the winds to inject on fully refined cells. - JW
                         ! Note this has to reflect whether this cell if fully refined would
                         ! contain the particle compared to the ref_radius, since ref_radius is
                         ! at full refinement. - JW
                         if (rad/(2.**(gr_maxRefine-lrefine(b)-1)) .le. ref_radius) then !else if 
#ifdef debug
print*, "rad =", rad, dr_globalMe
print*, "ref =", ref_radius, dr_globalMe
#endif
                            p_found = .true.
                         end if
                      end do
                      
                   end do
                end do
             end do

             deallocate(xc)
             deallocate(yc)
             deallocate(zc)

          end if      ! nodetype

          if (p_found) then
#ifdef debug
 print*, "[Particles_markRefineDerefine]: Particle found!", dr_globalMe
#endif
             if (lrefine(b) .lt. gr_maxRefine) then
                refine   (b) = .TRUE.
                derefine (b) = .FALSE.
                stay     (b) = .TRUE.
             endif
             if (lrefine(b) .eq. gr_maxRefine) then
                derefine (b) = .FALSE.
                stay     (b) = .TRUE.
             endif
          end if
       end do         ! loop over blocks
#endif
    deallocate(x, y, z, locx, locy, locz)
    return

end subroutine Particles_MarkRefineDerefine

!!****f* source/Particles/Particles_addNew
!!
!! NAME
!!    Particles_addNew
!!
!! SYNOPSIS
!!    call Particles_addNew( integer(in)  :: count,
!!                  optional,real(in)     :: pos(MDIM,count),
!!                           logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    This routine allows particles to be added during evolution.
!!    In the particles data structure it always initializes the tag and
!!    processor ID. If the optional argument "pos" is present then it
!!    will also initialize the position and block ID attributes in the
!!    particles. It returns the value FALSE in success if there isn't
!!    enough space left in the data structure for the requested number
!!    of particles.
!!
!! ARGUMENTS
!!
!!     count   :: the count of particles to be added
!!     pos     :: optional, contains the coordinates of the particles
!!     success :: This arg returns TRUE if there was enough space 
!!                in the particles data structure to add the requested
!!                number of particles, FALSE otherwise.
!!
!!    
!!  NOTES
!!
!! The constant MDIM is defined in constants.h .
!!
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_addNew (count, x, y, z, num_part_made, local_tags, ptype_in, mass_in, create_time)
  
  use Particles_Data, only : particles, pt_numLocal, pt_globalMe, &
            pt_maxPerProc, pt_meshComm, pt_meshMe, pt_indexList, &
            pt_indexCount, number_new_massive, new_massive_tags
  use Particles_sinkData, only : local_tag_number ! Global tags for both particles and particles_local - JW
  use Grid_interface, only : Grid_getBlkIDFromPos
  use Driver_data, only : dr_simTime
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

#define get_tag(arg1,arg2) ((arg1)*65536 + (arg2))
#define get_pno(arg1) ((arg1)/65536)
#define get_ppe(arg1) ((arg1) - get_pno(arg1)*65536)

  integer, INTENT(in)                          :: count
  real, dimension(count), intent(IN)           :: x, y, z
  real, optional, dimension(count), intent(IN) :: mass_in, create_time
  real, optional, intent(IN)                   :: ptype_in
  integer, intent(OUT)                         :: num_part_made
  real, dimension(count), intent(OUT)          :: local_tags


  real                   :: pt, ptype
  real, dimension(count) :: mass, c_time
  integer                :: block_no, part_PE, ire, i, pno

  logical, parameter     :: debug = .false.

if (present(mass_in)) then
    mass = mass_in
else
    mass = 0.0
end if

if (present(create_time)) then
    c_time = create_time
else
    c_time = dr_simTime
end if

#ifdef TYPE_PART_PROP
if (present(ptype_in)) then
    ptype = ptype_in
else
    ptype = 1.0
end if
#endif

pno = 0
num_part_made = 0

!print*, "Pos of ", 1, "star=", x(1), y(1), z(1), "on", pt_globalMe

do i=1, count

! Now get the proper processor and block for this location. If we are
! on the correct proc, we'll make the particle here.

  !print*, "Pos of ", i, "star=", x(i), y(i), z(i), "on", pt_globalMe

  call Grid_getBlkIDFromPos([x(i), y(i), z(i)], block_no, part_PE, pt_meshComm)

  if (part_PE .eq. pt_globalMe) then

    pt_numLocal = pt_numLocal + 1
    if (pt_numLocal .GT. pt_maxPerProc) &
       call Driver_abortFlash('Particles_addNew: Particle number exceeds pt_maxPerProc. Increase.')

    num_part_made = num_part_made + 1
    pno = pt_numLocal

    particles(:, pno) = 0.

    local_tag_number = local_tag_number + 1

! Note: If we form it using this procs local_tag_number, we can't cause
!       conflicts involving tags even if this particle belongs on a different
!       processor / block.
    particles(TAG_PART_PROP,pno)  = get_tag(local_tag_number, part_PE)
  

    particles(PROC_PART_PROP,pno)          = part_PE
    particles(BLK_PART_PROP,pno)           = block_no
    particles(POSX_PART_PROP,pno)          = x(i)
    particles(POSY_PART_PROP,pno)          = y(i)
    particles(POSZ_PART_PROP,pno)          = z(i)
    particles(CREATION_TIME_PART_PROP,pno) = c_time(i)
    particles(VELX_PART_PROP,pno)          = 0.0
    particles(VELY_PART_PROP,pno)          = 0.0
    particles(VELZ_PART_PROP,pno)          = 0.0
    particles(X_ANG_PART_PROP,pno)         = 0.0
    particles(Y_ANG_PART_PROP,pno)         = 0.0
    particles(Z_ANG_PART_PROP,pno)         = 0.0
    particles(MASS_PART_PROP,pno)          = mass(i)
    particles(ACCR_RATE_PART_PROP,pno)     = 0.0
    particles(ACCX_PART_PROP,pno)          = 0.0
    particles(ACCY_PART_PROP,pno)          = 0.0
    particles(ACCZ_PART_PROP,pno)          = 0.0
    particles(OACX_PART_PROP,pno)          = 0.0
    particles(OACY_PART_PROP,pno)          = 0.0
    particles(OACZ_PART_PROP,pno)          = 0.0
    particles(DTOLD_PART_PROP,pno)         = 0.0
#ifdef TYPE_PART_PROP
    particles(TYPE_PART_PROP,pno)          = ptype
#endif

    local_tags(i) = particles(TAG_PART_PROP, pno)

  ! Additions for AMUSE to learn if new particles were created during
  ! the Flash evolution step. Note we keep track of this variable in
  ! Particles_sinkData -JW
  
    number_new_massive = number_new_massive + 1
    new_massive_tags(number_new_massive) = int(particles(TAG_PART_PROP, pno))
  
  !print*, "Number of new sinks =", number_new_sinks
  !print*, "New tags = ", new_tags(1:number_new_sinks)
  
!  else
!    print*, 'Particles_addNew: WARNING! Particle loc proc /= to formation proc.'
!    call Driver_abortFlash('Particles_addNew: Particle loc proc /= to formation proc.')
  end if
  
end do

! Can't do this unless you call Particles_addNew from every proc.

! Now move all the particles to the proper processors and blocks.
!call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,&
!     pt_numLocal,pt_indexList,pt_indexCount,.true.)

return
  
end subroutine Particles_addNew

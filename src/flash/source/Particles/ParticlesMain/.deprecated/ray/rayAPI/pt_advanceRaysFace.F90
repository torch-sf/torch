!!
!! NAME
!!
!!  Particles_rayAdvance
!!
!! SYNOPSIS
!!
!!  Particles_rayAdvance(real(in)	:: dt)
!!
!! ARGUMENTS
!!  
!!   dt 						: timestep
!!  
!! DESCRIPTION
!!
!! main driver of raytracing algorithm
!! changes the particles data structure
!! also does MPI communication
!! loop over all particles and trace/split rays
!! assumes that all unsplit rays are adjacent in a block in the memory, also that memory is used
!! tightly so new rays are just appended to the end of the memory and processed in the
!! following loop
!!
!! COMMENTS
!!
!!
!!***

!#define DEBUG_RAY
!#define DEBUG_SANITY
#define PROC 0

! RAYS HAVE TO HAVE LAST TYPE ID
subroutine pt_advanceRaysFace(dt)

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
#include "GridParticles.h"

	include "Flash_mpi.h"

  real, intent(in) :: dt

	return
end subroutine pt_advanceRaysFace

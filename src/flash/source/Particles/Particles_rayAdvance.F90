!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014
!!see ray paper (ref. here)

!! Description:
!!   wrapper for raytracing routines
!!   calls generation and advancing routines
!!
!! Input: 
!!   dt: current simulation timestep

!! NOTES:
!!   ALL HARDCODED FOR 3D
!!   PERIODIC BOUNDARY CONDITIONS ARE FINE FOR DEFAULT SETTINGS, IF 
!!   FOR SOME REASON YOU WANT RAYS TO CROSS MULTIPLE TIMES LOOK IN pt_reconstructRay.f90
!!   Ghostbusters I > Ghostbusters II
!!
!! TODO figure out how to treat source in fully molecular medium
!!			I-front is generated one step after source turns on, as no atomic Hydrogen exists (DONE)
!! 
!! TODO disentangle rays from rest of the particle data structure, no rays survive after this, i.e.
!!      no stopped rays, i.e. c is infinity, otherwise persistent allocation needed (DONE)
!!
!! TODO more efficient MPI: instead of one unsorted array for all neighbouring processors,
!!      create 6 for all neighbouring ones, and populate those until some critical number of rays
!!      is filled, then do one point to point communication -> async (DONE)
!!
!! TODO rewrite MPI routine to be non-blocking in send and receive (DONE)
!!
!! TODO write some byte packing routine for MPI transport (DONE)
!!	
!! TODO split ray array into real and int part locally, keep MPI transport buffer as real, to not 
!!      setup two calls for int an real data (DONE)
!!
!! TODO switch out array data structure to doubly linked list 
!!
!! TODO internal threading each photon is independent, careful about creating new photons
!!      multiple lists that are joined together later? Seems like a better way than GPU
!!
!! TODO in the long run accessing in pt_solve is more efficient, for more detailed physics
!!
!!
!! -------------------------------------------------------------------
!! NON essential stuff
!!
!! TODO figure out what to copy to possible GPU extension, which fields. Ray generation on GPU or 
!!      locally?
!!
!! TODO explore alternating chemistry/ raytracing step for multiple energy bins:
!!      Lyman Werner -> Chemistry -> ionising rad. -> Chemistry -> next step 
!!
!! TODO check if reallocation each time step takes up time, alternatively switch to static array
!!      this increases memory footpring outside the raytracing routine though 
!!
!! TODO Disentangle communication routines from rest of the routines properly
!!
!! TODO switch out ray array with ray doubly linked list
!!
!! TODO at the moment only one type of rays (point or face) are able to be advanced
!!			depending on the loaded module the other source type is a stub
!!
!! TODO uses sources generated from WindDriving at the moment, decouple by writing some interfaces
!!      
!! play this when in use, crank it
!! Green Velvet & Harvard Bass - Lazer Beams

subroutine Particles_rayAdvance(dt)

#include "Flash.h"

  implicit none

  real, intent(in) :: dt

  return
end subroutine Particles_rayAdvance

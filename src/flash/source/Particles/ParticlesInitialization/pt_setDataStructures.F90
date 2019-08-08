!!****if* source/Particles/ParticlesInitialization/pt_setDataStructures
!!
!! NAME
!!
!!  pt_setDataStructures
!!
!! SYNOPSIS
!!
!!  call pt_setDataStructures()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

#include "Flash.h"
#include "Particles.h"
#include "constants.h"


subroutine pt_setDataStructures()
  use Particles_data, ONLY :  pt_typeInfo, pt_posInitialized, pt_velInitialized, pt_posAttrib,&
       pt_posPredAttrib, pt_velNumAttrib, pt_velAttrib, pt_velPredAttrib, pt_meshMe, particles, &
       therm_numAttrib, therm_Attrib

  use Simulation_interface, ONLY : Simulation_mapStrToInt,Simulation_mapParticlesVar
  use pt_interface, ONLY : pt_mapStringParamToInt, pt_picInit
  implicit none


  pt_posInitialized=.false.
  pt_velInitialized=.false.

  !! these data structures are introduced to make the interfaces uniform for the routines
  !!  where these attributes are to be used with different dimensionalities. With these
  !!  the different dimensionality can be managed through loops rather than ifdefs.
  pt_posAttrib(IAXIS)=POSX_PART_PROP
  pt_posAttrib(JAXIS)=POSY_PART_PROP
  pt_posAttrib(KAXIS)=POSZ_PART_PROP

  pt_posPredAttrib(:) = 1
#ifdef POSPREDX_PART_PROP
  pt_posPredAttrib(IAXIS)=POSPREDX_PART_PROP
  pt_posPredAttrib(JAXIS)=POSPREDY_PART_PROP
  pt_posPredAttrib(KAXIS)=POSPREDZ_PART_PROP
#endif

  pt_velNumAttrib=0
  pt_velAttrib=0
#ifdef VELX_PART_PROP
#ifdef VELX_VAR
  pt_velAttrib(PART_DS_IND,1)=VELX_PART_PROP
  pt_velAttrib(GRID_DS_IND,1)=VELX_VAR
  pt_velNumAttrib = 1

#ifdef VELY_PART_PROP
#ifdef VELY_VAR
  if(NDIM>1) then
     pt_velAttrib(PART_DS_IND,2)=VELY_PART_PROP
     pt_velAttrib(GRID_DS_IND,2)=VELY_VAR
     pt_velNumAttrib = 2
  else
     if (allocated(particles)) particles(VELY_PART_PROP,:)=0.0
  end if
#endif
#endif
  
#ifdef VELZ_PART_PROP
#ifdef VELZ_VAR
  if(NDIM>2) then
     pt_velAttrib(PART_DS_IND,3)=VELZ_PART_PROP
     pt_velAttrib(GRID_DS_IND,3)=VELZ_VAR
     pt_velNumAttrib = NDIM
  else
     if (allocated(particles)) particles(VELZ_PART_PROP,:)=0.0
  end if
#endif
#endif
#else
#ifdef DEBUG_PARTICLES
  if(pt_meshMe==MASTER_PE)print*,"no velocity attributes defined for particles"
#endif
#endif
#endif

#ifdef DEBUG_PARTICLES
  if(pt_meshMe==MASTER_PE)print*,"Particles_init: pt_velNumAttrib is", pt_velNumAttrib
  if(pt_meshMe==MASTER_PE)print*,"Particles_init: pt_velAttrib is", pt_velAttrib
#endif

pt_velPredAttrib(:,:)=1
#ifdef VELPREDX_PART_PROP
#ifdef VELX_VAR
  pt_velPredAttrib(PART_DS_IND,1)=VELPREDX_PART_PROP
  pt_velPredAttrib(GRID_DS_IND,1)=VELX_VAR
  
#ifdef VELPREDY_PART_PROP
#ifdef VELY_VAR
  if(NDIM>1) then
     pt_velPredAttrib(PART_DS_IND,2)=VELPREDY_PART_PROP
     pt_velPredAttrib(GRID_DS_IND,2)=VELY_VAR
  else
     if (allocated(particles)) particles(VELPREDY_PART_PROP,:)=0.0
  end if
#endif
#endif
  
#ifdef VELPREDZ_PART_PROP
#ifdef VELZ_VAR
  if(NDIM>2) then
     pt_velPredAttrib(PART_DS_IND,3)=VELPREDZ_PART_PROP
     pt_velPredAttrib(GRID_DS_IND,3)=VELZ_VAR
  else
     if (allocated(particles)) particles(VELPREDZ_PART_PROP,:)=0.0
  end if
#endif
#endif
#ifdef DEBUG_PARTICLES
  if(pt_meshMe==MASTER_PE)print*,"Particles_init: pt_velPredAttrib is", pt_velPredAttrib
#endif
#endif
#endif

! to track thermodynamical evolution in sb particle module uses overloaded particle properties

#ifdef TSN_PART_PROP
#ifdef MASS_PART_PROP
  therm_numAttrib = 2
  therm_Attrib(PART_DS_IND,1)= NSN_PART_PROP
  therm_Attrib(PART_DS_IND,2)= MASS_PART_PROP
  therm_Attrib(GRID_DS_IND,1)= TEMP_VAR
  therm_Attrib(GRID_DS_IND,2)= DENS_VAR
#endif
#endif

end subroutine pt_setDataStructures

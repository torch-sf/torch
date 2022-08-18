!!****if* source/Particles/ParticlesForces/longRange/gravity/ParticleMesh/Particles_longRangeForce
!!
!! NAME
!!
!!  Particles_longRangeForce
!!
!! SYNOPSIS
!!
!!  Particles_longRangeForce(real,intent(inout) :: particles(:,:),
!!                           integer,intent(in) :: p_count,
!!                           integer,intent(in) :: mapType
!!                           integer, intent(in), optional :: pot_var)
!!
!! DESCRIPTION
!!
!!  Computes long-range forces on particles, ie. forces which couple all
!!  particles to each other.  This version is for particle-mesh gravitation and
!!  maps the gravitational acceleration on the mesh to the particle positions.
!!  
!! ARGUMENTS
!!
!!   particles :: the list of particles to be operated on
!!   p_count   :: count of the particles in the list
!!   mapType   :: when mapping grid quantities to particle, method to use
!!   pot_var   :: if present, use this potential variable instead of GPOT_VAR
!!
!! NOTES
!!
!!  If GPOT_VAR and GPOT_PART_PROP are both defined, this subroutine
!!  additionally maps the gpot UNK variable to the gpot particle property.
!!
!!  Added optional potential variable argument. - JW 2017
!!  
!!***

!=======================================================================

subroutine Particles_longRangeForce (particles,p_count,mapType, pot_var)
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_mapMeshToParticles
  use Gravity_interface, ONLY : Gravity_accelListOfBlocks
  use Grid_interface, ONLY : Grid_mapMeshToParticles
  use Particles_data, ONLY : pt_posAttrib,pt_typeInfo

!-------------------------------------------------------------------------------

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
!-------------------------------------------------------------------------
  integer, intent(IN) :: p_count,mapType
  integer, intent(IN), optional :: pot_var
  real,dimension(NPART_PROPS,p_count),intent(INOUT) :: particles
  integer :: blockCount
  integer,parameter :: part_props=NPART_PROPS
  integer :: numAttrib
  integer :: potl
  integer,dimension(2,1) :: attrib

  integer,dimension(MAXBLOCKS) :: blockList

  if (.not. present(pot_var)) then
    potl = GPOT_VAR
  else
    potl = pot_var
  end if

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
! Map gravitational acceleration to particle positions

  numAttrib=1
  particles(ACCX_PART_PROP,:) = 0.
  particles(ACCY_PART_PROP,:) = 0.
  particles(ACCZ_PART_PROP,:) = 0.

  attrib(GRID_DS_IND,numAttrib)=GRAC_VAR
  attrib(PART_DS_IND,numAttrib)=ACCX_PART_PROP
  call Gravity_accelListOfBlocks(blockCount, blockList,IAXIS,GRAC_VAR, potl)

  call Grid_mapMeshToParticles(particles,part_props, BLK_PART_PROP,&
       p_count,pt_posAttrib,numAttrib,attrib,mapType)

  if (NDIM >= 2) then
     call Gravity_accelListOfBlocks(blockCount,blockList,JAXIS,GRAC_VAR, potl)
     attrib(PART_DS_IND,numAttrib)=ACCY_PART_PROP
     call Grid_mapMeshToParticles(particles,part_props, BLK_PART_PROP,&
          p_count,pt_posAttrib,numAttrib,attrib,mapType)
  endif

  if (NDIM == 3) then
     call Gravity_accelListOfBlocks(blockCount,blockList,KAXIS,GRAC_VAR, potl)
     attrib(PART_DS_IND,numAttrib)=ACCZ_PART_PROP
     call Grid_mapMeshToParticles(particles,part_props, BLK_PART_PROP,&
          p_count,pt_posAttrib,numAttrib,attrib,mapType)
  endif

#ifdef GPOT_VAR
#ifdef GPOT_PART_PROP
if (.not. present(pot_var)) then
  attrib(GRID_DS_IND,numAttrib)=GPOT_VAR
  attrib(PART_DS_IND,numAttrib)=GPOT_PART_PROP
  call Grid_mapMeshToParticles(particles,part_props, BLK_PART_PROP,&
       p_count,pt_posAttrib,numAttrib,attrib,mapType)
end if
#endif
#endif
 
!---------------------------------------------------------------------

  return

end subroutine Particles_longRangeForce

!=====================================================================

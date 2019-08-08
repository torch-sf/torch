!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_getFluxData
!!
!! NAME
!!  Grid_getFluxData
!!
!! SYNOPSIS
!!
!!
!!  call Grid_getFluxData(integer(IN) :: blockID,
!!                   integer(IN) :: axis,
!!                   real(INOUT) :: fluxes(NFLUXES,dataSize(1),dataSize(2),dataSize(3)),
!!                   integer(IN) :: dataSize(3),
!!          OPTIONAL,integer(IN) :: pressureSlots,
!!          OPTIONAL,real(IN)    :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!
!!  Get the fluxes in a direction specified by axis for boundary cells
!!  for block blockID. This routine needs to be used when using adaptive mesh
!!  since fluxes calculated by the two blocks that are at a fine/coarse
!!  boundary have different accuracy.
!!  
!!  This should be called after Grid_conserveFluxes, which gets the 
!!  fluxes updated for cells at fine/coarse boundaries and makes them consistent.
!!
!! ARGUMENTS
!!
!!  blockID : The local blockid
!!
!!
!!  axis : integer value specifying on which cell faces to get fluxes. 
!!         The options are IAXIS, JAXIS, or KAXIS defined in constants.h
!!
!!
!!  fluxes :  real array with space for fluxes, through one axis, 
!!            for all cells of a block and for all flux variables.
!!            fluxes(VAR, i, j, k) is VAR's flux through 
!!            the left cell face for cell i, j, k.  The 
!!            fluxes of the boundary cells of coarse blocks
!!            will have been appropriately changed, if the flux matching
!!            between fine and coarse boundaries (normally by a call to
!!            Grid_conserveFluxes), has occured since the last call
!!            of Grid_putFluxData.
!!
!!
!!  dataSize : integer array specifying the dimensions for fluxes
!!
!!             dataSize (1) holds the number of cells returned in the i direction
!!
!!             dataSize (2) holds the number of cells returned in the j direction
!!                          if 1 d problem, set datasize(2) = 1
!!
!!             dataSize (3) holds the number of cells returned in the k direction
!!                          if 1 or 2 d problem, set datasize(3) = 1
!!
!!             fluxes should contain space for fluxes of all cells in the block, 
!!             including guardcells.
!!
!!  pressureSlots: an optional array of integer values indicating fluxes that are to be
!!                 handled like pressures, rather than genuine fluxes.  If present, each
!!                 element greater than zero indicates one flux variable in the fluxes
!!                 array that may need special handling because it really does not scale
!!                 like a flux. Normally this would be pressure, or (e.g., with multiTemp
!!                 physics) several partial pressures given at cell interfaces. The
!!                 pressureSlots array could also list other flux variables that the
!!                 caller keeps in flux density form.
!!
!!  areaLeft :     areas of left and right faces, only used if special scaling is
!!                 requested with the pressureSlots argument.
!!
!! NOTES 
!!
!!   Any code calling this subroutine needs to know the explicit interface,
!!   since this interface contains optional dummy arguments and assumed-shape
!!   dummy arrays. Calling FORTRAN units should therefore contain a line like
!!       use Grid_interface, ONLY: Grid_getFluxData
!!
!!   This implementation is specific to Paramesh 4.
!!
!! SEE ALSO
!!
!!   Grid_putFluxData
!!   Grid_conserveFluxes
!!   hy_sweep
!!***

!! MODIFIED VERSION FOR THE MC TRACER PARTICLES
!! WRITES CORRECTED BOUNDARY FLUXES TO ITS RESPECTIVE FIELDS
!! BOUNDARY FLUX CORRECTION PERFORMED IN Grid_conserveFlux
!! ASSUMES GET FLUXES IS ONLY CALLED ONCE
!! IT WOULD BE BETTER TO APPLY BOUNDARY CORRECTION DIRECTLY 
!! DURING CALCULATION BUT THE FUNCTION IS PART OF PARAMESH


!!REORDER(5): flux_[xyz], gr_[xyz]flx
!!REORDER(4): fluxes
!!REORDER(5): gr_xflx_[yz]face, gr_yflx_[xz]face, gr_zflx_[xy]face
#include "Flash.h"
subroutine Grid_getFluxData(blockID, axis, fluxes, dataSize, pressureSlots, areaLeft)

  use physicaldata, ONLY : flux_x, flux_y, flux_z, nfluxes
  use tree, ONLY : surr_blks, nodetype
  use Grid_data, ONLY : gr_xflx, gr_yflx, gr_zflx

  use Grid_interface, ONLY : Grid_getBlkPtr,    &
                             Grid_releaseBlkPtr

#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >=2
  use Grid_data, ONLY : gr_xflx_yface, gr_yflx_xface
#if NDIM == 3
  use Grid_data, ONLY : gr_xflx_zface, gr_yflx_zface, gr_zflx_xface, gr_zflx_yface
#endif
#endif
#endif
  implicit none

!#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: blockID
  integer, intent(IN) :: axis
  integer, intent(IN), dimension(3) :: dataSize
  real, intent(INOUT), dimension(NFLUXES,dataSize(1),dataSize(2),dataSize(3)) :: fluxes
  integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
  real, intent(IN), OPTIONAL :: areaLeft(:,:,:)

#if NFLUXES > 0
  integer :: presVar, np
  integer,save,dimension(1),target :: presDefault = (/-1/)
  integer,pointer,dimension(:) :: presP
  integer :: sx,ex,sy,ey,sz,ez
  real, pointer, dimension(:,:,:,:) :: blockData

  if (present(pressureSlots)) then
     presP => pressureSlots
  else
     presP => presDefault
  end if

  ! Get block pointer
  call Grid_getBlkPtr(blockID,blockData,CENTER)

  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
  ex = dataSize(1)-NGUARD
  ey = dataSize(2)-NGUARD*K2D
  ez = dataSize(3)-NGUARD*K3D

  select case(axis)
  case(IAXIS)

#ifndef CHOMBO_COMPATIBLE_HYDRO
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
     fluxes(:,sx:ex+1,sy:sy+1,sz:ez) = gr_xflx_yface(:,:,1:2,:,blockID)
     fluxes(:,sx:ex+1,ey:ey+1,sz:ez) = gr_xflx_yface(:,:,3:4,:,blockID)
#if NDIM == 3
     fluxes(:,sx:ex+1,sy:ey,sz:sz+1) = gr_xflx_zface(:,:,:,1:2,blockID)
     fluxes(:,sx:ex+1,sy:ey,ez:ez+1) = gr_xflx_zface(:,:,:,3:4,blockID)
#endif
#endif
#endif
#endif
     fluxes(:,sx,  sy:ey,sz:ez) = flux_x(:nfluxes,1,:,:,blockID)
     fluxes(:,ex+1,sy:ey,sz:ez) = flux_x(:nfluxes,2,:,:,blockID)
     fluxes(:,sx+1,sy:ey,sz:ez) = gr_xflx(:,1,:,:,blockID) 
     fluxes(:,ex,  sy:ey,sz:ez) = gr_xflx(:,2,:,:,blockID)

#ifdef DENS_FLUX
 		 blockData(XMFLX_VAR,sx,	sy:ey,sz:ez) = flux_x(DENS_FLUX,1,:,:,blockID)
		 blockData(XMFLX_VAR,ex+1,sy:ey,sz:ez) = flux_x(DENS_FLUX,2,:,:,blockID)
 		 blockData(XMFLX_VAR,sx+1,sy:ey,sz:ez) = gr_xflx(DENS_FLUX,1,:,:,blockID) 
 		 blockData(XMFLX_VAR,ex,  sy:ey,sz:ez) = gr_xflx(DENS_FLUX,2,:,:,blockID)
#endif

!     print*,sx,sy,sz
!     print*,ex,ey,ez
!     print*,fluxes(DENS_VAR,:,sy,sz)
!     print*,flux_x(:,1:2,sy,sz,blockID),blockID
!     print*,gr_xflx(:,1:2,sy,sz,blockID),blockID

     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,1,1+K2D,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,1,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx,sy:ey,sz:ez).NE.0.0) &
                   fluxes(presVar,sx,sy:ey,sz:ez) = fluxes(presVar,sx,sy:ey,sz:ez) / areaLeft(sx,sy:ey,sz:ez)
           end if
           if (.NOT.(surr_blks(1,3,1+K2D,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,3,1+K2D,1+K3D,blockID) == nodetype(blockID))) then
              fluxes(presVar,ex+1,sy:ey,sz:ez) = fluxes(presVar,ex+1,sy:ey,sz:ez) / areaLeft(ex+1,sy:ey,sz:ez)
           end if
        end if
     end do
#ifdef CHOMBO_COMPATIBLE_HYDRO
     fluxes(:,sx,  sy:ey,sz:ez) = fluxes(:,sx,  sy:ey,sz:ez) - fluxes(:,sx+1,sy:ey,sz:ez) !new - old
     fluxes(:,sx+1,sy:ey,sz:ez) = 0.0
     fluxes(:,ex+1,sy:ey,sz:ez) = fluxes(:,ex+1,sy:ey,sz:ez) - fluxes(:,ex,  sy:ey,sz:ez) !new - old
     fluxes(:,ex,  sy:ey,sz:ez) = 0.0
#endif

  case(JAXIS)
#ifndef CHOMBO_COMPATIBLE_HYDRO
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM >= 2
     fluxes(:,sx:sx+1,sy:ey+1,sz:ez) = gr_yflx_xface(:,1:2,:,:,blockID)
     fluxes(:,ex:ex+1,sy:ey+1,sz:ez) = gr_yflx_xface(:,3:4,:,:,blockID)
#if NDIM == 3
     fluxes(:,sx:ex,sy:ey+1,sz:sz+1) = gr_yflx_zface(:,:,:,1:2,blockID)
     fluxes(:,sx:ex,sy:ey+1,ez:ez+1) = gr_yflx_zface(:,:,:,3:4,blockID)
#endif
#endif
#endif
#endif
     fluxes(:,sx:ex,sy,  sz:ez) = flux_y(:nfluxes,:,1,:,blockID) 
     fluxes(:,sx:ex,ey+1,sz:ez) = flux_y(:nfluxes,:,2,:,blockID) 
     fluxes(:,sx:ex,sy+1,sz:ez) = gr_yflx(:,:,1,:,blockID) 
     fluxes(:,sx:ex,ey,  sz:ez) = gr_yflx(:,:,2,:,blockID)

#ifdef DENS_FLUX
 		 blockData(YMFLX_VAR,sx:ex,sy,  sz:ez) = flux_y(DENS_FLUX,:,1,:,blockID) 
		 blockData(YMFLX_VAR,sx:ex,ey+1,sz:ez) = flux_y(DENS_FLUX,:,2,:,blockID)
 		 blockData(YMFLX_VAR,sx:ex,sy+1,sz:ez) = gr_yflx(DENS_FLUX,:,1,:,blockID) 
 		 blockData(YMFLX_VAR,sx:ex,ey,  sz:ez) = gr_yflx(DENS_FLUX,:,2,:,blockID)
#endif

#if NDIM > 1
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
           if (.NOT.(surr_blks(1,2,1,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,1,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx:ex,sy,sz:ez).NE.0.0) &
                fluxes(presVar,sx:ex,sy,sz:ez) = fluxes(presVar,sx:ex,sy,sz:ez) / areaLeft(sx:ex,sy,sz:ez)
           end if
           if (.NOT.(surr_blks(1,2,3,1+K3D,blockID) > 0 .AND. &
             surr_blks(3,2,3,1+K3D,blockID) == nodetype(blockID))) then
              where (areaLeft(sx:ex,ey+1,sz:ez).NE.0.0) &
                   fluxes(presVar,sx:ex,ey+1,sz:ez) = fluxes(presVar,sx:ex,ey+1,sz:ez) / areaLeft(sx:ex,ey+1,sz:ez)
           end if
        end if
     end do
#endif
#ifdef CHOMBO_COMPATIBLE_HYDRO
     fluxes(:,sx:ex,sy,  sz:ez) = fluxes(:,sx:ex,sy,  sz:ez) - fluxes(:,sx:ex,sy+1,sz:ez)
     fluxes(:,sx:ex,sy+1,sz:ez) = 0.0
     fluxes(:,sx:ex,ey+1,sz:ez) = fluxes(:,sx:ex,ey+1,sz:ez) - fluxes(:,sx:ex,ey,  sz:ez)
     fluxes(:,sx:ex,ey,  sz:ez) = 0.0
#endif

  case(KAXIS)
#ifndef CHOMBO_COMPATIBLE_HYDRO
#ifdef FLASH_HYDRO_UNSPLIT
#if NDIM == 3
     fluxes(:,sx:sx+1,sy:ey,sz:ez+1) = gr_zflx_xface(:,1:2,:,:,blockID)
     fluxes(:,ex:ex+1,sy:ey,sz:ez+1) = gr_zflx_xface(:,3:4,:,:,blockID)
     fluxes(:,sx:ex,sy:sy+1,sz:ez+1) = gr_zflx_yface(:,:,1:2,:,blockID)
     fluxes(:,sx:ex,ey:ey+1,sz:ez+1) = gr_zflx_yface(:,:,3:4,:,blockID)
#endif
#endif
#endif
     fluxes(:,sx:ex,sy:ey,sz  ) = flux_z(:nfluxes,:,:,1,blockID)
     fluxes(:,sx:ex,sy:ey,ez+1) = flux_z(:nfluxes,:,:,2,blockID)
     fluxes(:,sx:ex,sy:ey,sz+1) = gr_zflx(:,:,:,1,blockID)
     fluxes(:,sx:ex,sy:ey,ez  ) = gr_zflx(:,:,:,2,blockID)

#ifdef DENS_FLUX
 		 blockData(ZMFLX_VAR,sx:ex,sy:ey,sz  ) = flux_z(DENS_FLUX,:,:,1,blockID)
		 blockData(ZMFLX_VAR,sx:ex,sy:ey,ez+1) = flux_z(DENS_FLUX,:,:,2,blockID)
 		 blockData(ZMFLX_VAR,sx:ex,sy:ey,sz+1) = gr_zflx(DENS_FLUX,:,:,1,blockID)
 		 blockData(ZMFLX_VAR,sx:ex,sy:ey,ez  ) = gr_zflx(DENS_FLUX,:,:,2,blockID)
#endif

#if NDIM > 2
     do np = 1,size(presP,1)
        presVar = presP(np)
        if (presVar > 0) then
!!$        where (areaLeft(sx:ex,sy:ey,sz).NE.0.0) &  ! should not happen in any supported geometry
           if (.NOT.(surr_blks(1,2,2,1,blockID) > 0 .AND. &
             surr_blks(3,2,2,1,blockID) == nodetype(blockID))) then
              fluxes(presVar,sx:ex,sy:ey,sz) = fluxes(presVar,sx:ex,sy:ey,sz) / areaLeft(sx:ex,sy:ey,sz)
           end if
           if (.NOT.(surr_blks(1,2,2,3,blockID) > 0 .AND. &
             surr_blks(3,2,2,3,blockID) == nodetype(blockID))) then
              fluxes(presVar,sx:ex,sy:ey,ez+1) = fluxes(presVar,sx:ex,sy:ey,ez+1) / areaLeft(sx:ex,sy:ey,ez+1)
           end if
        end if
     end do
#endif
#ifdef CHOMBO_COMPATIBLE_HYDRO
     fluxes(:,sx:ex,sy:ey,sz) = fluxes(:,sx:ex,sy:ey,sz) - fluxes(:,sx:ex,sy:ey,sz+1)
     fluxes(:,sx:ex,sy:ey,sz+1) = 0.0
     fluxes(:,sx:ex,sy:ey,ez+1) = fluxes(:,sx:ex,sy:ey,ez+1) - fluxes(:,sx:ex,sy:ey,ez)
     fluxes(:,sx:ex,sy:ey,ez) = 0.0
#endif

  end select
#endif


  call Grid_releaseBlkPtr(blockID,blockData,CENTER)

  return
end subroutine Grid_getFluxData

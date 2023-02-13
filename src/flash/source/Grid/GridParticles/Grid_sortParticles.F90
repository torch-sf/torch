!!****if* source/Grid/GridParticles/Grid_sortParticles
!!
!! NAME
!!  Grid_sortParticles
!!
!! SYNOPSIS
!!
!!  call Grid_sortParticles(real(INOUT)    :: dataBuf(:,:),
!!                          integer(IN)    :: props
!!                          integer(INOUT) :: localCount,
!!                          integer(IN)    :: elementTypes,
!!                          integer(IN)    :: maxCount,
!!                          integer(OUT)   :: elementsPerBlk(:,:),
!!                          integer(IN)    :: attrib1,
!!                 optional,integer(IN)    :: attrib2)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Sorts the elements by block number. There are three types of blocknumber associated with
!!  elements data structures which have valid values, but are not valid blocknumber in the 
!!  mesh. These are "UNKNOWN", "LOST" and "NONEXISTENT". The sorter finds out the number of blocks 
!!  on the current processor in the mesh, and puts all elements associated with blocknumbers 
!!  in the range of 1 to localNumBlocks in the processor into the appropriate bins. In the bin
!!  for localNumBlocks+1 it puts all elements with block number = UNKNOWN, and in 
!!  localNumBloks+2 it put all elements with block number = LOST, in localNumBlocks+3
!!  it puts elements with block number NONEXISTENT. For any other block 
!!  number in the "gr_ptBlk" field of any element, the routine aborts. 
!!
!!  If block number pointer is provided in attrib1 and particle types pointer in attrib2 but
!!  the particles have not been assigned a type, this routine will crash giving a
!!  free(): invalid size error. - SCL 
!!
!! ARGUMENTS 
!!
!!  dataBuf : List of elements. It is two dimensional real array, the first dimension
!!              represents each element's properties, and second dimension is index to
!!              elements.
!!
!!  props : number of properties of each element in the dataBuf datastructure
!!
!!  localCount : While coming in it contains the current number of elements mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local elements might change, and the new value is put back into it.
!!  elementTypes  : Count of different types of elements in the simulation
!!  maxCount : This is parameter determined at runtime, and is the maximum number of local
!!                        elements that a simulation expects to have. All the arrays that hold
!!                        particles in the Particles unit are allocated based on this number.
!!                        
!! elementsPerBlk : integer array containing number of elements on each blk.
!!  attrib1       : the primary property of the elements on which to sort them
!!  attrib2       : if present, then the elements are first sorted on attrib2, and then 
!!                  within each group with a given attrib2 value, on the attrib1.
!!
!!
!! NOTES
!!   Currently this routine is called (only?) by io_writeParticles in permanent guardcell mode,
!!   by Particles_advance, by Particles_updateGridVar, and by the routines that map mesh to particles
!!   or particles to mesh.
!!   CAUTION : This routine should not be called when " gr_ptSourceBuf" needs to have valid values
!!             since it is used as scratch space by this routine.
!!
!!   The algorithm used for sorting requires (and silently assumes) that
!!     o   if attrib1 or attrib2 is BLK_PART_PROP then for each element it is 
!!         either in the valid range 1..MAXBLOCKS,
!!         or has one of the special values NONEXISTENT, UNKNOWN, LOST; and
!!     o   if the attribute in particles is the particle type
!!         property for each element, it is in the valid range 1..elementTypes.
!!   Elements with a block value of UNKNOWN or LOST will be effectively dropped (by
!!   potentially reusing their storage locations for other elements, and by not counting them
!!   in the 'elementsPerBlk' output array), but this routine does *not* update the localCount
!!   counter; so the caller better take care of updating its view of the number of elements (probably
!!   from the information in 'elementsPerBlk') after Grid_sortElements returns, if some of the
!!   elements may be NONEXISTENT.
!!   Elements with a block value of NONEXISTENT will be effectively dropped (by
!!   potentially reusing their storage locations for other elements, and by not counting them
!!   in the 'elementsPerBlk' output array), and this routine *does* update the localCount
!!   counter by subtracting the number of such elements.
!!
!! DEV: This seems only correctly implemented for attrib1=BLK_PART_PROP, attrib2=TYPE_PART_PROP(or not present) 
!!
!! VorAMR + Torch Note:
!! If block number pointer is provided in attrib1 and particle types pointer in attrib2 but
!! the particles have not been assigned a type, this routine will crash giving a
!! free(): invalid size error. - SCL
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_SORT_PARTICLES
#endif

subroutine Grid_sortParticles(dataBuf,props, localCount,elementTypes,maxCount,&
                              elementsPerBlk,attrib1, attrib2)
  use gr_ptData, ONLY : gr_ptSourceBuf, gr_ptLogLevel, gr_ptBlk, gr_ptKeepLostParticles
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Particles.h"
  integer,intent(IN) ::  maxCount,props,elementTypes
  integer, intent(INOUT) :: localCount
  integer, parameter :: UNKNOWNBLK=1, LOSTBLK=2, NONEXISTENTBLK=3, &
       ALLBLOCKS=MAXBLOCKS+3
  real,dimension(props,maxCount),intent(INOUT) :: dataBuf
  integer,dimension(MAXBLOCKS,elementTypes),intent(OUT) :: elementsPerBlk
  integer,dimension(ALLBLOCKS,elementTypes)::localElementsPerBlk
  integer, intent(IN) :: attrib1
  integer,optional, intent(IN) :: attrib2
  integer, dimension((MAXBLOCKS)*elementTypes+3) :: pntr
  integer :: i, j,k,m,n,attrib,localNumBlocks,validParticles,lastBlkID
  integer :: grptBlkProp
  logical :: pType
  
  integer, dimension(elementTypes) :: unknownBin, lostBin, nonexistentBin
  integer :: validPointer, unknownID, lostID, nonexistentID

  localElementsPerBlk(:,:)=0

  validPointer=0
#ifdef DEBUG_SORT_PARTICLES
  print*,'Grid_sortParticles local/maxCount=',localCount,maxCount,',attrib1=',attrib1,'attrib2=',attrib2
#endif
  
  pType = present(attrib2)
  attrib = 1
  if(pType) attrib = attrib2 

  !! If no sorting is to be done on type, and there is only one block
  !! then all elements belong to the same block and no more sorting
  !! needs to be done

!  if((.not.pType).and.(localNumBlocks==1)) then
!     localElementsPerBlk(1,1)=localCount
!  else

     !! either there is more than one type, or more than one block, so sort
     !!  This is not an inplace sort, and so is order N in complexity
     !!  There are two passes, one to find out the number of elements
     !!  in each block, and the second one to put elements in order
     
     !! First copy all elements to a scrach, and count number 
     !! of elements per block  
     
     
  k=1 
  call Grid_getLocalNumBlks(localNumBlocks)
  unknownID=localNumBlocks+UNKNOWNBLK
  lostID=localNumBlocks+LOSTBLK
  nonexistentID=localNumBlocks+NONEXISTENTBLK
    
  do i = 1,localCount
     if (pType) k=int(dataBuf(attrib,i)) !DEV: Should protect against attrib > elementTypes here! - KW
     j = int(dataBuf(attrib1,i))
#ifdef DEBUG_PARTICLES
     if((j>localNumBlocks).or.(j<LOST)) then
        call Driver_abortFlash("Grid_sortParticles : undefined block number")
     end if
#endif
     validPointer=validPointer+1
     gr_ptSourceBuf(1:props,validPointer)=dataBuf(1:props,i)
     if(j==UNKNOWN)j=unknownID
     if(j==LOST)j=lostID
     if(j==NONEXISTENT)j=nonexistentID
     localElementsPerBlk(j,k)=localElementsPerBlk(j,k)+1
  end do
  
  localCount=validPointer
  !! Setup starting point for each block as it would
  !! once they are sorted
  
  pntr(1)=1
  k=2

  do j=1,elementTypes
     do i = 1,localNumBlocks
        pntr(k)=pntr(k-1)+localElementsPerBlk(i,j)
        k=k+1
     end do
  end do
  
  do i = localNumBlocks+1,nonexistentID-1
     pntr(k)=pntr(k-1)+sum(localElementsPerBlk(i,1:elementTypes))
     k=k+1
  end do
  
  !! Drop elements in their right place
  k=0

  do i = 1,localCount
     j=int(gr_ptSourceBuf(attrib1,i))

     if(j>0) then
        if(pType)k=int(gr_ptSourceBuf(attrib,i))-1 !DEV: Should protect against attrib > elementTypes here! - KW
        n=k*localNumBlocks+j
     else
        n=localNumBlocks*elementTypes
        if(j==UNKNOWN)n=n+UNKNOWNBLK
        if(j==LOST)n=n+LOSTBLK
        if(j==NONEXISTENT)n=n+NONEXISTENTBLK
     end if
     m=pntr(n)
     dataBuf(:,m)=gr_ptSourceBuf(:,i)
     pntr(n)=pntr(n)+1
  end do

  elementsPerBlk(1:MAXBLOCKS,:)=localElementsPerBlk(1:MAXBLOCKS,:)
  localCount=localCount-sum(localElementsPerBlk(nonexistentID,1:elementTypes))

  return
end subroutine Grid_sortParticles

!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_sanitizeDataAfterInterp
!!
!! NAME
!!
!!  gr_sanitizeDataAfterInterp
!!
!!
!! SYNOPSIS
!!
!!  call gr_sanitizeDataAfterInterp(integer(in)          :: blkList(count),
!!                                  integer(in)          :: count,
!!                                  character(len=*)(in) :: info,
!!                                  integer(in)          :: layers(MDIM))
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  check whether solution data in certain variables lie in a reasonable
!!  range of values.
!!
!!  Energies (ENER_VAR and EINT_VAR) are expected to be .ge. gr_smalle,
!!  and the density (DENS_VAR) is expected to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!!  For data that do satisfy these expectations, warning messages are
!!  generated, but the offending data is not modified.
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   layers - number of guardcell layers to be included in the check 
!!
!! NOTES
!!
!!  The checks are based on gr_conserveToPrimitive, which is called to
!!  convert solution data back from conserved form when using the old
!!  (convertToConsvdForMeshCalls) way of ensuring that the mesh
!!  handles data interpolation in conserved form.
!!
!!  This is meant to be called where gr_conserveToPrimitive used to be
!!  called when using the new (convertToConsvdInMeshInterp) way of
!!  ensuring that the mesh handles data interpolation in conserved
!!  form.
!!
!! SEE ALSO
!!
!!  gr_conserveToPrimitive
!!
!!
!! BUGS
!!
!!  This routine accesses the global variable storage 
!!  array unk directly.  It won't work for data stored
!!  in the paramesh workspace array WORK. It won't work
!!  for the Uniform Grid (its functionality is currently
!!  not needed there). 
!!
!!***

!!REORDER(5): unk

#define DEBUG_CONSCONV

subroutine gr_sanitizeDataAfterInterp(blkList,count, info, layers)

  use Grid_data, ONLY : gr_smallrho,gr_smalle, gr_meshMe, &
       gr_sanitizeDataMode, gr_sanitizeVerbosity
  use Logfile_interface, ONLY : Logfile_stamp
  use physicaldata, ONLY:unk, gcell_on_cc
  use tree, ONLY:nodetype, lrefine
  use paramesh_dimensions, ONLY: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd, kl_bndi, ndim

  implicit none
#include "constants.h"
#undef REAL_FORMAT
#define REAL_FORMAT "(1PG23.16)"
#include "Flash.h"

  integer,intent(IN) :: count
  integer, dimension(count), intent(IN) :: blkList
  character(len=*), intent(IN) :: info
  integer,dimension(MDIM), intent(IN):: layers
  integer :: n, block

  integer ::  i,j
  integer :: iskip, jskip, kskip
  integer :: il,iu,jl,ju,kl,ku
  integer :: kwrite,locs(3),kReorder(1:ku_bnd-kl_bnd+1),nReorder
  character(len=32), dimension(5,2) :: block_buff
  character(len=32)                 :: number_to_str

111 format (a,a,a1,(1x,a18,'=',a),(1x,a2,'=',a5),(1x,a5,'=',a),(1x,a4,'=',a))
112 format (i3,1x,24(1x,1G8.2))
113 format (' :,',i2,',',i2,1x,24(1x,1G8.2))

  if (gr_sanitizeDataMode == 0) return ! IMMEDIATE RETURN

  iskip = NGUARD - layers(IAXIS)
  jskip = (NGUARD - layers(JAXIS)) * K2D
  kskip = (NGUARD - layers(KAXIS)) * K3D
  il = il_bnd + iskip
  iu = iu_bnd - iskip
  jl = jl_bnd + jskip
  ju = ju_bnd - jskip
  kl = kl_bnd + kskip
  ku = ku_bnd - kskip

  nReorder = 0

  do n = 1,count
     block=blkList(n)

#ifdef DENS_VAR
     if (gcell_on_cc(DENS_VAR) .OR. gr_sanitizeDataMode == 2) then
        ! small limits -- in case the interpolants are not monotonic
        if (any(unk(DENS_VAR,il:iu,jl:ju,kl:ku,block) .LT. gr_smallrho)) then
           kwrite = kl_bndi
           if (gr_sanitizeVerbosity .GE. 5 .AND. ndim==3) then
              call set_kReorder
!              print*,'kReorder(1:nReorder)',kReorder(1:nReorder)
              locs = minloc(unk(DENS_VAR,il:iu,jl:ju,kReorder(1:nReorder),block))
!              print*,'LOCS:',locs
              kwrite = kReorder(locs(3))
           end if
           write (block_buff(1,1), '(a18)') 'min. unk(DENS_VAR)'
           !        write (number_to_str, '('//REAL_FORMAT//',a1)') minval(unk(DENS_VAR,il:iu,jl:ju,kl:ku,block)), ','
           write (number_to_str, '(G30.22)') minval(unk(DENS_VAR,il:iu,jl:ju,kl:ku,block))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') gr_meshMe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') block
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'type'
           write (number_to_str, '(i7)') nodetype(block)
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

		   write (block_buff(5,1), '(a)') 'ref.'
		   write (number_to_str, '(i7)') lrefine(block)
		   write (block_buff(5,2), '(a)') trim(adjustl(number_to_str))

           if (gr_sanitizeVerbosity .GE. 1) call Logfile_stamp( block_buff, 4, 2, 'WARNING '//info)
           if (gr_sanitizeVerbosity .GE. 4) print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,5)

           if (gr_sanitizeVerbosity .GE. 5) then
              do j=ju_bnd,jl_bnd,-1
                 if (kwrite==kl_bndi) then
                 ! For 3D, this prints a slice at the lowest k index that is interior - KW
                    print 112, j, (unk(DENS_VAR,i,j,kwrite,block), i=il_bnd,iu_bnd)
                 else
                    print 113, j,kwrite, (unk(DENS_VAR,i,j,kwrite,block), i=il_bnd,iu_bnd)
                 end if
              end do
           end if

           if (gr_sanitizeDataMode == 3) then
              unk(DENS_VAR,il:iu,jl:ju,kl:ku,block) = max(gr_smallrho,unk(DENS_VAR,il:iu,jl:ju,kl:ku,block))
           end if
           if (gr_sanitizeDataMode == 4) call Driver_abortFlash("DENS var below acceptable minimum")
        end if
     end if
     
#endif

#ifdef ENER_VAR               
     if (gcell_on_cc(ENER_VAR) .OR. gr_sanitizeDataMode == 2) then
        ! energy
        if (any(unk(ENER_VAR,il:iu,jl:ju,kl:ku,block) .LT. gr_smalle*0.999999999)) then
           kwrite = kl_bndi
           if (gr_sanitizeVerbosity .GE. 5 .AND. ndim==3) then
              call set_kReorder
              locs = minloc(unk(ENER_VAR,il:iu,jl:ju,kReorder(1:nReorder),block))
              kwrite = kReorder(locs(3))
           end if
           write (block_buff(1,1), '(a)') 'min. unk(ENER_VAR)'
           write (number_to_str, '('//REAL_FORMAT//')') minval(unk(ENER_VAR,il:iu,jl:ju,kl:ku,block))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') gr_meshMe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') block
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'type'
           write (number_to_str, '(i7)') nodetype(block)
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

		   write (block_buff(5,1), '(a)') 'ref.'
		   write (number_to_str, '(i7)') lrefine(block)
		   write (block_buff(5,2), '(a)') trim(adjustl(number_to_str))

           if (gr_sanitizeVerbosity .GE. 1) call Logfile_stamp( block_buff, 4, 2, 'WARNING '//info)
           if (gr_sanitizeVerbosity .GE. 4) print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,5)

           if (gr_sanitizeVerbosity .GE. 5) then
              do j=ju_bnd,jl_bnd,-1
                 if (kwrite==kl_bndi) then
                    print 112, j, (unk(ENER_VAR,i,j,kwrite,block), i=il_bnd,iu_bnd) 
                 else
                    print 113, j,kwrite, (unk(ENER_VAR,i,j,kwrite,block), i=il_bnd,iu_bnd)
                 end if
              end do
           end if

           if (gr_sanitizeDataMode == 3) then
              unk(ENER_VAR,il:iu,jl:ju,kl:ku,block) = max(gr_smalle,unk(ENER_VAR,il:iu,jl:ju,kl:ku,block))
           end if
           if (gr_sanitizeDataMode == 4) call Driver_abortFlash("ENER var below acceptable minimum")
        end if
     end if
#endif
#ifdef EINT_VAR
     if (gcell_on_cc(EINT_VAR) .OR. gr_sanitizeDataMode == 2) then
        if (any(unk(EINT_VAR,il:iu,jl:ju,kl:ku,block) .LT. gr_smalle*0.999999999)) then
           kwrite = kl_bndi
           if (gr_sanitizeVerbosity .GE. 5 .AND. ndim==3) then
              call set_kReorder
              locs = minloc(unk(EINT_VAR,il:iu,jl:ju,kReorder(1:nReorder),block))
              kwrite = kReorder(locs(3))
           end if
           write (block_buff(1,1), '(a)') 'min. unk(EINT_VAR)'
           write (number_to_str, '('//REAL_FORMAT//')') minval(unk(EINT_VAR,il:iu,jl:ju,kl:ku,block))
           write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(2,1), '(a)') 'PE'
           write (number_to_str, '(i7)') gr_meshMe
           write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(3,1), '(a)') 'block'
           write (number_to_str, '(i7)') block
           write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

           write (block_buff(4,1), '(a)') 'type'
           write (number_to_str, '(i7)') nodetype(block)
           write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

		   write (block_buff(5,1), '(a)') 'ref.'
		   write (number_to_str, '(i7)') lrefine(block)
		   write (block_buff(5,2), '(a)') trim(adjustl(number_to_str))

           if (gr_sanitizeVerbosity .GE. 1) call Logfile_stamp( block_buff, 4, 2, 'WARNING '//info)
           if (gr_sanitizeVerbosity .GE. 4) print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,5)

           if (gr_sanitizeVerbosity .GE. 5) then
              do j=ju_bnd,jl_bnd,-1
                 if (kwrite==kl_bndi) then
                    print 112, j, (unk(EINT_VAR,i,j,kwrite,block), i=il_bnd,iu_bnd) 
                 else
                    print 113, j,kwrite, (unk(EINT_VAR,i,j,kwrite,block), i=il_bnd,iu_bnd)
                 end if
              end do
           end if

           if (gr_sanitizeDataMode == 3) then
              unk(EINT_VAR,il:iu,jl:ju,kl:ku,block) = max(gr_smalle,unk(EINT_VAR,il:iu,jl:ju,kl:ku,block))
           end if
           if (gr_sanitizeDataMode == 4) call Driver_abortFlash("EINT var below acceptable minimum")
        end if
     end if
#endif

  end do

  return 

contains
  subroutine set_kReorder
    integer :: i,j,k
    if (nReorder==0) then       !We need to do this only once per gr_sanitize... call.

       !  The following code sets
       !  kReorder = (/5,6,7,8,9,10,11,12,4,13,3,14,2,15,1,16/) ! for kl:ku = 1:16

       i = 0
       do k = kl, ku
          if (k>NGUARD .AND. k .LE. ku_bnd-NGUARD) then
             i = i+1
             kReorder(i) = k
          end if
       end do
       do j = 1, layers(KAXIS)
          k = NGUARD+1-j
          if (k .LE. NGUARD) then
             i = i+1
             kReorder(i) = k
          end if
          k = ku_bnd-NGUARD+j
          if (k > ku_bnd-NGUARD) then
             i = i+1
             kReorder(i) = k
          end if
       end do
       nReorder = i
    end if
  end subroutine set_kReorder

end subroutine gr_sanitizeDataAfterInterp
        

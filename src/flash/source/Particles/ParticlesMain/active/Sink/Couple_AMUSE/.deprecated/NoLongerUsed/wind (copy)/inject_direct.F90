
subroutine inject_direct(loc, injectMass, injectVelocity, twind, dt)

!#define debug

use Grid_data, ONLY: gr_meshComm, gr_meshMe

use Grid_interface, ONLY: Grid_getBlkIDFromPos, Grid_getBlkPhysicalSize, &
    Grid_getBlkNeighBlkIDFromPos, Grid_getDeltas, Grid_getCellCoords, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
    Grid_getBlkCenterCoords, Grid_fillGuardCells, Grid_getMinCellSize, &
    Grid_notifySolnDataUpdate
    
use Driver_interface, ONLY : Driver_abortFlash
    
use Eos_interface, ONLY : Eos_wrapped

use pt_windInterface, only : overlap

#include "Flash.h"
#include "constants.h"

implicit none

#include "Flash_mpi.h"

integer, parameter :: dp = kind(1.d0)

real(dp), intent(in)    :: loc(3)
real(dp), intent(in)    :: injectMass, injectVelocity, twind
real(dp), intent(inout) :: dt

real(dp) :: star_x, star_y, star_z, cell_pos(3)
real(dp) :: local_pos(3), blockCenter(3), blockSize(3), delta(3)
real(dp) :: cell_ipos, cell_jpos, cell_kpos, cell_top(3), cell_bot(3)
real(dp) :: dVol


real(dp) :: injectThermE, injectRadius, injectDens, injectKinE, injectTotE
real(dp) :: overlap_frac,  sumMass
real(dp) :: idir, jdir, kdir, rad, dx, dy, dz, xvel, yvel, zvel

real(dp) :: background_density, t_0, R_1, kinFrac

real(dp), allocatable, dimension(:) :: xc, yc, zc
integer, allocatable, dimension(:) :: affectedBlks

integer :: size_x, size_y, size_z
integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM), blockID
integer :: i, j, k, ii, jj, kk, ii1, jj1, kk1, ilocal, jlocal, klocal, n
integer :: Proc_ID, proclocal, blocklocal, ind_array(3)
integer :: n_cell_depth, numBlks
integer :: neighProc, neighBlk, neighInd(3), neighDir(3)

integer, parameter :: trans_tag = 103180
integer, parameter :: mess_tag  = 801031
integer :: messages, trans_status(MPI_STATUS_SIZE), ierr

real(dp), pointer, dimension(:,:,:,:) :: solndata

real(dp), parameter :: yr = 60_dp*60_dp*24_dp*365_dp

logical :: off_proc, off_block, iHaveInjectBlk, dupBlk, inside_cell

!!!! Debugging

real(dp) :: before(3), after(3)

!!!!

n_cell_depth = 2
iHaveInjectBlk = .false.
numBlks=0
sumMass = 0.0_dp

call Grid_getCellIndexFromPosition(Proc_ID, blockID, &
                               ii, jj, kk, &
                               loc(1), loc(2), loc(3))
                               

ind_array = [ ii, jj, kk ]

! Get average background density in the injection area.

background_density = 0.0_dp

if (gr_meshMe == Proc_ID) then

  call Grid_getBlkPtr(blockID, solndata)

  do kk1=kk-n_cell_depth, kk+n_cell_depth
    do jj1=jj-n_cell_depth, jj+n_cell_depth
      do ii1=ii-n_cell_depth, ii+n_cell_depth

        background_density = background_density  & 
                           + solndata(DENS_VAR, ii1, jj1, kk1)
                           
      end do
    end do
  end do
  
  call Grid_releaseBlkPtr(blockID, solndata)
end if

background_density = background_density / ((2.0_dp*real(n_cell_depth)+1.0_dp)**3.0_dp)

call MPI_BCAST(ind_array , 3, MPI_INT, Proc_ID, gr_meshComm, ierr)
call MPI_ALLREDUCE(MPI_IN_PLACE, background_density, 1, MPI_REAL, MPI_SUM, gr_meshComm, ierr)

write(*,'(A,ES13.3e3)') "background density=", background_density
                               
ii = ind_array(1); jj = ind_array(2); kk = ind_array(3)

                              

!call Grid_getDeltas(blockID, delta)
call Grid_getMinCellSize(delta(1))
delta=delta(1)
dVol = delta(1)*delta(2)*delta(3)

injectRadius = sqrt(3.0_dp)*delta(1)

!n_cell_depth = 2
iHaveInjectBlk = .false.
numBlks=0
sumMass = 0.0_dp

injectThermE = 0.0_dp

! Figure out where the edge of the free streaming cold kinetic wind
! is by calculating where the mass of the wind = mass of swept up material.

! Having second thoughts about this idea, since we want to actually see
! this evolve because of the kinetic expansion. I think what I actually
! want to know is does the wind blown inside the cell with the star
! ever not appear 100% kinetic? If it does, do we later evolve to a
! time when it *is* 100% kinetic?

t_0 = sqrt(3.0_dp * (injectMass / dt) / (4.0_dp * PI * background_density * injectVelocity**3.0_dp))

R_1 = injectVelocity * t_0

write(*, '(A,ES13.3e3)') "Free expansion edge R_1 from simple mass calculation = ", R_1

! Weaver et al 77 eqn 12
R_1 = .74296_dp*(1.0_dp / background_density * injectMass / dt)**(3.0_dp/10.0_dp) &
    * injectVelocity**(1.0_dp/10.0_dp) * twind**(2.0_dp/5.0_dp)

write(*, '(A,ES13.3e3)') "Free expansion edge R_1 from Weaver 77 = ", R_1
write(*, '(A,ES13.3e3)') "t_0 =", t_0
write(*, '(A,ES13.3e3)') "twind =", twind

! Calculate the injection density assuming that the injection mass
! is spread equally around a sphere of radius=injection Radius.

injectDens = injectMass/(4.0_dp/3.0_dp*PI*injectRadius**3.0_dp)

! Calculate the total amount of energy that would be injected into
! a cell if 100% was kinetic with the injectDens and velocity
! for the minimum cell size.

injectTotE = 0.5_dp * injectDens * dVol * injectVelocity**2.0_dp

!! Get the physical locations of all the cells on this block. Necessary for
!! getting the neighbor block location below.

call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)

size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

allocate(xc(size_x))
allocate(yc(size_y))
allocate(zc(size_z))

! n_cell_depth + 1 for the central cell.
allocate(affectedBlks((n_cell_depth+1)**3))

do kk1=kk-n_cell_depth, kk+n_cell_depth
  do jj1=jj-n_cell_depth, jj+n_cell_depth
    do ii1=ii-n_cell_depth, ii+n_cell_depth
      
      local_pos(1) = loc(1) + (ii1 - ii)*delta(1)
      local_pos(2) = loc(2) + (jj1 - jj)*delta(2)
      local_pos(3) = loc(3) + (kk1 - kk)*delta(3)
      
      call Grid_getCellIndexFromPosition(Proc_ID, blockID, &
                                       ilocal, jlocal, klocal, &
                                       local_pos(1), local_pos(2), local_pos(3))
                                       
      if (gr_meshMe .eq. Proc_ID) then                                
                                       
        call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xc, size_x)
        call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yc, size_y)
        call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zc, size_z)
        
        cell_top = [xc(ilocal), yc(jlocal), zc(klocal)] + 0.5_dp*delta
        cell_bot = [xc(ilocal), yc(jlocal), zc(klocal)] - 0.5_dp*delta
                                       
    ! Calculate the fractional directions from the center of the sphere
    ! to the cell center (for the velocity directions).
    
        dx = xc(ilocal) - loc(1)
        dy = yc(jlocal) - loc(2)
        dz = zc(klocal) - loc(3)
    
        rad = sqrt(dx**2.0_dp + dy**2.0_dp + dz**2.0_dp)
        
        !injectDens = injectMass /(4.0_dp * PI * rad**2.0_dp)
      
        if (rad .ne. 0.0_dp) then
    
          idir = dx/rad
          jdir = dy/rad
          kdir = dz/rad
        
        else
      
          idir = 0.0_dp
          jdir = 0.0_dp
          kdir = 0.0_dp
      
        endif

        
        ! Calculate the fractional position within the cell of the
        ! wind source between the cell wall and the cell center in
        ! each direction.
        
        inside_cell = .false.
        
        cell_ipos = abs(dx/(0.5_dp*delta(1)))
        cell_jpos = abs(dy/(0.5_dp*delta(2)))
        cell_kpos = abs(dz/(0.5_dp*delta(3)))
        
        ! If we are between the midpoint of the cell and the cell wall
        ! the momentum injection should be scaled appropriately in each
        ! direction. So if we are at the wall, momentum should be 100%.
        ! At the cell center, it should be 0%.
        
        if (cell_ipos .lt. 1.0_dp) idir = cell_ipos*idir
                    
        if (cell_jpos .lt. 1.0_dp) jdir = cell_jpos*jdir
        
        if (cell_kpos .lt. 1.0_dp) kdir = cell_kpos*kdir
        
        if ((cell_ipos .lt. 1.0_dp) .and. &
            (cell_jpos .lt. 1.0_dp) .and. &
            (cell_kpos .lt. 1.0_dp)) inside_cell = .true.
        
        ! Calculate the velocity in each direction, given the direction
        ! of the radius vector and whether the source is inside this cell.
        
        xvel = idir * injectVelocity
        yvel = jdir * injectVelocity
        zvel = kdir * injectVelocity
        
        ! Now calculate the overlapping areas of the sphere and this
        ! cell.
                                   
        call overlap(1, injectRadius, loc, cell_bot, &
                  cell_top, 100, overlap_frac)

                  
        if (overlap_frac .gt. 0.0_dp) then

#ifdef debug
        ! Print normed velocity.
        write(*,'(A,ES13.3e3)') "Normed velocity = ", sqrt(xvel**2.0_dp + yvel**2.0_dp + zvel**2.0_dp)
        if (abs(sqrt(xvel**2.0_dp + yvel**2.0_dp + zvel**2.0_dp) - injectVelocity) .gt. 1d-6) then
          write(*,'(A,ES13.3e3)') " inject velocity = ", injectVelocity
        stop
        end if
#endif
        ! Note R_1 is the wind shock edge as defined in Weaver et al 1977.
          if (rad .ge. R_1) then
            ! If your past the wind shock, 
            ! divide the energy equally b/t thermal and
            ! kinetic (its roughly correct, Weaver et al 77 and Koo and McKee 92).
            kinFrac = 0.5_dp
            ! Also, the density is dominated by the swept up material
            ! in this region, so just inject that.
            ! Weaver et al 77 eqn 11.5.
            
            ! Maybe we just let the hydro code sort out where gas got swept up... 
            !injectDens = 1.1016_dp * (injectMass / dt) & 
            !           / ( PI * R_1**2.0 * injectVelocity)
                       
            ! Finally, velocity is the post shock velocity, not the injection velocity
            ! (Weaver eqn 11.4).
            xvel = idir * injectVelocity / 4.0_dp * 0.90773_dp * (R_1/rad)**2.0_dp
            yvel = jdir * injectVelocity / 4.0_dp * 0.90773_dp * (R_1/rad)**2.0_dp
            zvel = kdir * injectVelocity / 4.0_dp * 0.90773_dp * (R_1/rad)**2.0_dp
            
          else
            kinFrac = 1.0_dp
            print*, "Fully kinetic in here."
            
            ! Lets consider the mass equally spread around an expanding shell,
            ! such that the density goes like 1/r^2.
            ! Probably a better approximation than the density being constant
            ! everywhere.
            
            !injectDens = injectMass /(4.0_dp * PI * rad**2.0_dp)
            
          end if

          
          ! Calculate the kinetic energy actually injected into this cell.
          injectKinE = kinFrac * 0.5_dp * overlap_frac*injectDens*dVol &
                   * (xvel**2.0_dp + yvel**2.0_dp + zvel**2.0_dp)

          ! Is the source in this cell?
          
          if (inside_cell) then !(overlap_frac .eq. 1.0_dp) .and. 
            print*, "Source inside this cell, converting some &
                    kinetic energy into thermal energy."

            ! Calculate the difference between whats actually injected and
            ! what analytically should be injected. This we will inject as
            ! thermal energy instead.
          
            injectThermE = injectTotE - injectKinE
          
          else
          
            injectThermE = (1.0_dp - kinFrac)*injectKinE
            
          end if
        
          iHaveInjectBlk = .true.
          ! Store affected blocks.
          ! Check if this block is already stored.
          dupBlk = .false.
          do n=1, numBlks
            if(affectedBlks(n) == blockID) then
              dupBlk = .true.
              exit
            endif
          end do
          ! If not, store it now.
          if (dupBlk .eqv. .false.) then
            numBlks = numBlks + 1
            affectedBlks(numBlks) = blockID
          end if
      
          call Grid_getBlkPtr(blockID, solndata)
      
#ifdef debug          
          write(*,'(A, 3I4)') "ii, jj, kk", ii1, jj1, kk1
          write(*,'(A, 3I4)') "ilocal, jlocal, klocal", ilocal, jlocal, klocal 
          write(*,'(A,ES13.3E3)') "Overlap frac off blk =", overlap_frac
          write(*,'(A,ES13.3E3)') "rad off blk =", rad
          write(*,'(A,ES13.3E3)') "injectRadius off blk =", injectRadius
          write(*,'(A,3ES13.3E3)') "Cell bot =", cell_bot
          write(*,'(A,3ES13.3E3)') "Cell top =", cell_top
          write(*,'(A,3ES13.3E3)') "idir, jdir, kdir", idir, jdir, kdir
          write(*,*) "Kinetic energy =", injectKinE
          write(*,*) "Thermal energy =", injectThermE
          write(*,*) "Total energy =", injectTotE
          if (inside_cell) write(*,*) "Yo I'm in this cell!"

          before = solndata(VELX_VAR:VELZ_VAR, ilocal, jlocal, klocal)

          write(*,'(A,ES13.3E3)') "Solndata off blk dens before =", solndata(DENS_VAR, ilocal, jlocal, klocal)
          write(*,'(A,3ES13.3E3)') "Solndata off blk vel  before =", before
          call flush()
#endif

          solndata(VELX_VAR:VELZ_VAR, ilocal, jlocal, klocal) = &
           (solndata(DENS_VAR, ilocal, jlocal, klocal)*solndata(VELX_VAR:VELZ_VAR, ilocal, jlocal, klocal) & 
            + overlap_frac*injectDens*[xvel, yvel, zvel]) &
            / (solndata(DENS_VAR, ilocal, jlocal, klocal) + overlap_frac*injectDens)
          
          solndata(DENS_VAR, ilocal, jlocal, klocal) = & 
          solndata(DENS_VAR, ilocal, jlocal, klocal) + overlap_frac*injectDens
          
          solndata(EINT_VAR, ilocal, jlocal, klocal) = &
          solndata(EINT_VAR, ilocal, jlocal, klocal) + injectThermE/(overlap_frac*injectDens*dVol) ! Flash uses ergs / g
#ifdef debug
          after = solndata(VELX_VAR:VELZ_VAR, ilocal, jlocal, klocal)
          write(*,'(A,ES13.3E3)') "Solndata off blk dens after =", solndata(DENS_VAR, ilocal, jlocal, klocal)
          write(*,'(A,3ES13.3E3)') "Solndata off blk vel  after =", after
          call flush()
#endif
          call Grid_releaseBlkPtr(blockID, solndata)

!#ifdef debug
!          if (overlap_frac .ne. 0.0 .or. before(1) .ne. after(1)) stop
!#endif


!#ifdef debug
          sumMass = sumMass + overlap_frac*injectDens*dVol
!          print*, overlap_frac
!#endif
        end if
      end if
    end do
  end do
end do

#ifdef debug

if (abs(sumMass - injectMass)/sumMass .gt. 1.0d-4) then

  print*, "Oh noes! Total injectMass not actually injected!"
  print*, "injectMass =, sumMass ="
  print*,  injectMass, sumMass
  call Driver_abortFlash("[inject_direct]: Wind mass injected /= to wind mass from star!")

end if
#endif

deallocate(xc)
deallocate(yc)
deallocate(zc)

! Update the EOS variables if we have an affected block for only those blocks.

call Grid_notifySolnDataUpdate()

if (iHaveInjectBlk) then
  do n=1, numBlks
#ifdef debug
    print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc, numBlks ", affectedBlks(n), gr_meshMe, numBlks
#endif
    call Grid_getBlkIndexLimits(affectedBlks(n),blkLimits,blkLimitsGC)
    call Eos_wrapped(MODE_DENS_EI, blkLimits, affectedBlks(n))
  end do
!end if
else
#ifdef debug
   print*, "I'm proc", gr_meshMe, "and I have no injection blocks."
#endif
!  return
end if

deallocate(affectedBlks)

!call Grid_fillGuardCells(CENTER, ALLDIR)

!messages  = 0
!off_block = .false.
!off_proc  = .false.

!!loc = [star_x, star_y, star_z]

!call Grid_getBlkIDFromPos(loc, blockID, Proc_ID, gr_meshComm)

!#ifdef debug
!print*, "Comm = ", gr_meshComm
!print*, "MyProc = ", gr_meshMe
!print*, "Proc_ID = ", Proc_ID
!#endif

!if (gr_meshMe == Proc_ID) then

!  call Grid_getBlkCenterCoords(blockID, blockCenter)
!  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, CENTER)
!  call Grid_getBlkPhysicalSize(blockID, blockSize)
!  call Grid_getDeltas(blockID, delta)

!  injectRadius = 3.0*delta(1) !sqrt(3.0)*delta(1)

!!write(*,*) MDIM

!do i=1, MDIM
  
!  local_pos(i) = loc(i) - blockCenter(i) + blockSize(i)/2.0
  
!end do

!#ifdef debug
!  write(*,*) "loc=", loc 
!  write(*,*) "blockCenter=", blockCenter
!  write(*,*) 'blockSize=', blockSize
!  write(*,*) "local_pos=", local_pos
!  write(*,*) 'delta=', delta
!#endif
 
!      ii = floor(local_pos(1)/delta(1)) + 1 + NGUARD ! index in the block
!      cell_ipos = mod(local_pos(1),delta(1)) / delta(1) ! fractional loc in the cell

!  if (MDIM .gt. 1) then
!      jj = floor(local_pos(2)/delta(2)) + 1 + NGUARD
!      cell_jpos = mod(local_pos(2),delta(2)) / delta(2)
!  else
!      jj = 0
!      cell_jpos = 0.0
!  end if

!  if (MDIM .gt. 2) then
!      kk = floor(local_pos(3)/delta(3)) + 1 + NGUARD
!      cell_kpos = mod(local_pos(3),delta(3)) / delta(3)
!  else
!      kk = 0
!      cell_kpos = 0.0
!  end if

!#ifdef debug
!  write(*,*) "i,j,k =", ii,jj,kk
!  write(*,*) "local_pos =", local_pos
!  write(*,*) "block_id =", blockID
!#endif

!! Get the physical locations of all the cells on this block. Necessary for
!! getting the neighbor block location below.

!  size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
!  size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
!  size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

!  allocate(xc(size_x))
!  allocate(yc(size_y))
!  allocate(zc(size_z))

!  call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xc, size_x)
!  call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., yc, size_y)
!  call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., zc, size_z)

!  dVol = delta(1)*delta(2)*delta(3)

!! Okay now we know which cell the star is in (and where in the cell it is)
!! so lets drop some stuff in the grid. - JW

!! Number of cells to go in each direction. Note that a 3x3x3 injection
!! covers a 4x4x4 cube.

!  n_cell_depth = ceiling(injectRadius/delta(1))
!#ifdef debug
!  print*, "cell_depth =", n_cell_depth
!#endif
!  !stop

!do kk1=kk-n_cell_depth, kk+n_cell_depth
!  do jj1=jj-n_cell_depth, jj+n_cell_depth
!    do ii1=ii-n_cell_depth, ii+n_cell_depth

!! Have to be careful here... we are about
!! to put mass, vel etc into cells that
!! might be in guard cells or even on
!! another processor. - JW
!! So lets figure out which processor, block
!! and cell its on each loop. There's
!! probably a cheaper way, but its in
!! testing for now. - JW
      
!      off_block = .false.
!      off_proc  = .false.
!      before = 0.0
!      after  = 0.0
!      overlap_frac = 0.0
      
!      !x_cc = loc(1) + (ii1-ii)*delta(1)
!      !y_cc = loc(2) + (jj1-jj)*delta(2)
!      !z_cc = loc(3) + (kk1-kk)*delta(3)
      
!      if (ii1 .lt. blkLimits(LOW,IAXIS)) then
!        off_block  = .true. 
!        neighDir(1) = LEFT_EDGE
!        neighInd(1) = ii1 + NXB - NGUARD
!      else if (ii1 .gt. blkLimits(HIGH,IAXIS)) then
!        off_block = .true.
!        neighDir(1) = RIGHT_EDGE
!        neighInd(1) = ii1 - NXB + NGUARD
!      else
!        neighDir(1) = CENTER
!        neighInd(1) = ii1
!      end if
!      if (jj1 .lt. blkLimits(LOW,JAXIS)) then
!        off_block  = .true. 
!        neighDir(2) = LEFT_EDGE
!        neighInd(2) = jj1 + NYB - NGUARD
!      else if (jj1 .gt. blkLimits(HIGH,JAXIS)) then
!        off_block = .true.
!        neighDir(2) = RIGHT_EDGE
!        neighInd(2) = jj1 - NYB + NGUARD
!      else
!        neighDir(2) = CENTER
!        neighInd(2) = jj1
!      end if
!      if (kk1 .lt. blkLimits(LOW,KAXIS)) then
!        off_block  = .true. 
!        neighDir(3) = LEFT_EDGE
!        neighInd(3) = kk1 + NZB - NGUARD
!      else if (kk1 .gt. blkLimits(HIGH,KAXIS)) then
!        off_block = .true.
!        neighDir(3) = RIGHT_EDGE
!        neighInd(3) = kk1 - NZB + NGUARD
!      else
!        neighDir(3) = CENTER
!        neighInd(3) = kk1
!      end if

!      if (off_block) then
!!          print*, xc(ii1)
!!          print*, yc(jj1)
!!          print*, zc(kk1)
!!          print*, delta(1), delta(2), delta(3)
!!          print*, ii, jj, kk, ii1, jj1, kk1
!          cell_pos = [xc(ii1), yc(jj1), zc(kk1)]
!!          cell_pos = [loc(1)*(ii1*delta(1), loc(2)*jj1*delta(2), loc(3)*kk1*delta(3)]
!!          print*, cell_pos
!!          print*, "Hey, one of these cells is not on this block."
!          call Grid_getBlkNeighBlkIDFromPos(blockID, cell_pos, neighDir, neighBlk, neighProc)
!!          print*, "We are currently on proc ", Proc_ID, " and block ", blockID
!!          print*, "This cell is on proc ", neighProc, " and block ", neighBlk
!!          print*, "and its i, j, k are ", neighInd
          
!        if (Proc_ID .ne. neighProc) off_proc = .true.
      
!      end if
      
!      if ( off_block .and. off_proc) then
      
!        ! send info for this cell to the other processor
    
!        messages = 1
        
!        call MPI_Send(messages, 1, MPI_INTEGER, neighProc, mess_tag, &
!                    gr_meshComm, ierr)
      
!        call MPI_Send(neighBlk, 1, MPI_INTEGER, neighProc, trans_tag, &
!                    gr_meshComm, ierr)
                    
!        call MPI_Send(neighInd, 3, MPI_INTEGER, neighProc, trans_tag, &
!                    gr_meshComm, ierr)
        
!        ! and go on to the next loop for this processor
        
!        messages = 0
        
!        continue
      
!      else if (off_block .and. (.not. off_proc)) then
        
!        !deallocate(xc)
!        !deallocate(yc)
!        !deallocate(zc)
        
!        call Grid_getBlkIndexLimits(neighBlk, blkLimits, blkLimitsGC, CENTER)
!        call Grid_getDeltas(neighBlk, delta)
      
!        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
!        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
!        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

!        !allocate(xc(size_x))
!        !allocate(yc(size_y))
!        !allocate(zc(size_z))

!        call Grid_getCellCoords(IAXIS, neighBlk, CENTER, .true., xc, size_x)
!        call Grid_getCellCoords(JAXIS, neighBlk, CENTER, .true., yc, size_y)
!        call Grid_getCellCoords(KAXIS, neighBlk, CENTER, .true., zc, size_z)

!        dVol = delta(1)*delta(2)*delta(3)
        
!        ! Calculate the fractional directions from the center of the sphere
!        ! to the cell center (for the velocity directions).
        
!        dx = loc(1) - xc(neighInd(1))
!        dy = loc(2) - yc(neighInd(2))
!        dz = loc(3) - zc(neighInd(3))
        
!        rad = sqrt(dx**2.0 + dy**2.0 + dz**2.0)
        
!        idir = dx/rad
!        jdir = dy/rad
!        kdir = dz/rad
        
!        ! Get the fractional volume overlap of the sphere and the cell
        
!          call overlap(1, injectRadius, loc, [xc(neighInd(1))-0.5*delta(1), & 
!                       yc(neighInd(2))-0.5*delta(2), zc(neighInd(3))-0.5*delta(3)], &
!                      [xc(neighInd(1))+0.5*delta(1), yc(neighInd(2))+0.5*delta(2), & 
!                      zc(neighInd(3))+0.5*delta(3)], 100, overlap_frac)
!#ifdef debug        
!        write(*,'(A,ES13.3E3)') "Overlap frac =", overlap_frac
!        write(*,'(A,ES13.3E3)') "rad =", rad
!        write(*,'(A,ES13.3E3)') "injectRadius =", injectRadius
!        write(*,'(A,3ES13.3E3)') "Cell bot =", [xc(neighInd(1))-0.5*delta(1), & 
!                       yc(neighInd(2))-0.5*delta(2), zc(neighInd(3))-0.5*delta(3)]
!        write(*,'(A,3ES13.3E3)') "Cell top =",  [xc(neighInd(1))+0.5*delta(1), &
!                       yc(neighInd(2))+0.5*delta(2), zc(neighInd(3))+0.5*delta(3)]
!#endif

!        call Grid_getBlkPtr(neighBlk, solndata)
    
!#ifdef debug
!        before = solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3))
!        write(*,'(A,ES13.3E3)') "Solndata dens before =", solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))
!        write(*,'(A,3ES13.3E3)') "Solndata vel  before =", before
!        !write(*,'(A,ES13.3E3)') "Solndata eint before =", solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))
!#endif      
!        ! Inject velocity in a momentum conserving way.
!        solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3)) = &
!          (solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))*dVol*solndata(VELX_VAR:VELZ_VAR,neighInd(1),neighInd(2),neighInd(3)) &
!           + overlap_frac*injectMass*injectVelocity*[idir, jdir, kdir]) &
!           / (solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))*dVol + overlap_frac*injectMass)
!        ! Inject mass.
!        solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3)) = & 
!        solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3)) + overlap_frac*injectMass / dVol
!        ! Inject thermal energy.
!        solndata(EINT_VAR, neighInd(1),neighInd(2),neighInd(3)) = &
!        solndata(EINT_VAR, neighInd(1),neighInd(2),neighInd(3)) + overlap_frac*injectThermE/injectMass ! Flash uses ergs / g

!#ifdef debug        
!        write(*,'(A,ES13.3E3)') "Solndata dens after =", solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))
!        after = solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3))
!        write(*,'(A,3ES13.3E3)') "Solndata vel  after =", after
!        call flush()
!        if ((overlap_frac .ne. 0.0) .or. before(1) .ne. after(1)) stop
!#endif        
      
!        call Grid_releaseBlkPtr(neighBlk, solndata)
        
!        call Eos_wrapped(MODE_DENS_EI, blkLimits, neighBlk)
        
!      else
      
!        if (Proc_ID .eq. gr_meshMe) then
      
!          call Grid_getBlkPtr(blockID, solndata)
                  
!        ! Calculate the fractional directions from the center of the sphere
!        ! to the cell center (for the velocity directions).
        
!          dx = loc(1) - xc(ii1)
!          dy = loc(2) - yc(jj1)
!          dz = loc(3) - zc(kk1)
        
!          rad = sqrt(dx**2.0 + dy**2.0 + dz**2.0)
        
!          idir = dx/rad
!          jdir = dy/rad
!          kdir = dz/rad
        
!        ! Get the fractional volume overlap of the sphere and the cell
        
!          call overlap(1, injectRadius, loc, [xc(ii1)-0.5*delta(1), & 
!                       yc(jj1)-0.5*delta(2), zc(kk1)-0.5*delta(3)], &
!                      [xc(ii1)+0.5*delta(1), yc(jj1)+0.5*delta(2),  &
!                       zc(kk1)+0.5*delta(3)], 100, overlap_frac)
!#ifdef debug          
!          write(*,'(A, 3I4)') "ii, jj, kk", ii1, jj1, kk1             
!          write(*,'(A,ES13.3E3)') "Overlap frac off blk =", overlap_frac
!          write(*,'(A,ES13.3E3)') "rad off blk =", rad
!          write(*,'(A,ES13.3E3)') "injectRadius off blk =", injectRadius
!          write(*,'(A,3ES13.3E3)') "Cell bot =", [xc(ii1)-0.5*delta(1), & 
!                       yc(jj1)-0.5*delta(2), zc(kk1)-0.5*delta(3)]
!          write(*,'(A,3ES13.3E3)') "Cell top =",  [xc(ii1)+0.5*delta(1), &
!                                                  yc(jj1)+0.5*delta(2),  &
!                                                  zc(kk1)+0.5*delta(3)]

!          before = solndata(VELX_VAR:VELZ_VAR, ii1, jj1, kk1)

!          write(*,'(A,ES13.3E3)') "Solndata off blk dens before =", solndata(DENS_VAR, ii1, jj1, kk1)
!          write(*,'(A,3ES13.3E3)') "Solndata off blk vel  before =", before
!#endif

!          solndata(VELX_VAR:VELZ_VAR, ii1, jj1, kk1) = &
!           (solndata(DENS_VAR, ii1, jj1, kk1)*dVol*solndata(VELX_VAR:VELZ_VAR, ii1, jj1, kk1) & 
!            + overlap_frac*injectMass*injectVelocity*[idir, jdir, kdir]) &
!            / (solndata(DENS_VAR, ii1, jj1, kk1)*dVol + overlap_frac*injectMass)
          
!          solndata(DENS_VAR, ii1, jj1, kk1) = & 
!          solndata(DENS_VAR, ii1, jj1, kk1) + overlap_frac*injectMass / dVol
          
!          solndata(EINT_VAR, ii1, jj1, kk1) = &
!          solndata(EINT_VAR, ii1, jj1, kk1) + overlap_frac*injectThermE/injectMass ! Flash uses ergs / g
!#ifdef debug
!          after = solndata(VELX_VAR:VELZ_VAR, ii1, jj1, kk1)
!          write(*,'(A,ES13.3E3)') "Solndata off blk dens after =", solndata(DENS_VAR, ii1, jj1, kk1)
!          write(*,'(A,3ES13.3E3)') "Solndata off blk vel  after =", after
!          call flush()
!#endif
!          call Grid_releaseBlkPtr(blockID, solndata)

!#ifdef debug
!          if (overlap_frac .ne. 0.0 .or. before(1) .ne. after(1)) stop
!#endif

!          call Eos_wrapped(MODE_DENS_EI, blkLimits, blockID)
!#ifdef debug
!          if (xc(ii1)==delta(1) .and. yc(jj1)==delta(2) .and. zc(kk1)==delta(3)) then
!            write(*,*) "In central cell."
!            write(*,'(3ES13.3E3)') xc(ii1), yc(jj1), zc(kk1)
!            call flush()
!            stop
!          end if
!#endif
!        end if
        
!      end if

!    end do
!  end do
!end do

!! We're done here, so let all the procs know.

!deallocate(xc)
!deallocate(yc)
!deallocate(zc)

!messages = -1

!call MPI_Send(messages, 1, MPI_INTEGER, neighProc, mess_tag, &
!            gr_meshComm, ierr)


!else  ! I'm not the proc with the star on it, so listen for messages.

!  ! Wait for a call from the proc with the star telling you that
!  ! you are needed for something.

!  do while (messages .ge. 0)
  
!    call MPI_Recv(messages, 1, MPI_INTEGER, MPI_ANY_SOURCE, mess_tag, &
!             gr_meshComm, trans_status, ierr)
             
!    !print*, "Recieved message =", messages
             
!      if (messages .eq. 1) then
      
!        call MPI_Recv(neighBlk, 1, MPI_INTEGER, MPI_ANY_SOURCE, trans_tag, &
!             gr_meshComm, trans_status, ierr)
!        call MPI_Recv(neighInd, 3, MPI_INTEGER, MPI_ANY_SOURCE, trans_tag, &
!             gr_meshComm, trans_status, ierr)
        
!        call Grid_getBlkIndexLimits(neighBlk, blkLimits, blkLimitsGC, CENTER)
!        call Grid_getDeltas(neighBlk, delta)
        
!        injectRadius = sqrt(3.0)*delta(1)
      
!        size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
!        size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
!        size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

!        allocate(xc(size_x))
!        allocate(yc(size_y))
!        allocate(zc(size_z))

!        call Grid_getCellCoords(IAXIS, neighBlk, CENTER, .true., xc, size_x)
!        call Grid_getCellCoords(JAXIS, neighBlk, CENTER, .true., yc, size_y)
!        call Grid_getCellCoords(KAXIS, neighBlk, CENTER, .true., zc, size_z)

!        dVol = delta(1)*delta(2)*delta(3)
                
!        ! Calculate the fractional directions from the center of the sphere
!        ! to the cell center (for the velocity directions).
        
!        dx = loc(1) - xc(neighInd(1))
!        dy = loc(2) - yc(neighInd(2))
!        dz = loc(3) - zc(neighInd(3))
        
!        rad = sqrt(dx**2.0 + dy**2.0 + dz**2.0)
        
!        idir = dx/rad
!        jdir = dy/rad
!        kdir = dz/rad
        
!        ! Get the fractional volume overlap of the sphere and the cell
        
!        call overlap(1, injectRadius, loc, [xc(neighInd(1))-0.5*delta(1), &
!                     yc(neighInd(2))-0.5*delta(2), zc(neighInd(3))-0.5*delta(3)], &
!                    [xc(neighInd(1))+0.5*delta(1), yc(neighInd(2))+0.5*delta(2), &
!                     zc(neighInd(3))+0.5*delta(3)], 100, overlap_frac)
        
!        call Grid_getBlkPtr(neighBlk, solndata)
!#ifdef debug        
!        write(*,'(A,ES13.3E3)') "Solndata off proc dens before =", solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))
!        write(*,'(A,3ES13.3E3)') "Solndata off proc vel  before =", solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3))
!#endif
!        ! Inject velocity in a momentum conserving way.
!        solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3)) = &
!          (solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))*dVol*solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3)) &
!           + overlap_frac*injectMass*injectVelocity*[idir, jdir, kdir]) &
!           / (solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))*dVol + overlap_frac*injectMass)
!        ! Inject mass.
!        solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3)) = & 
!        solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3)) + overlap_frac*injectMass / dVol
!        ! Inject thermal energy.
!        solndata(EINT_VAR, neighInd(1),neighInd(2),neighInd(3)) = &
!        solndata(EINT_VAR, neighInd(1),neighInd(2),neighInd(3)) + overlap_frac*injectThermE/injectMass ! Flash uses ergs / g
!#ifdef debug
!        write(*,'(A,ES13.3E3)') "Solndata dens off proc after =", solndata(DENS_VAR, neighInd(1),neighInd(2),neighInd(3))
!        write(*,'(A,3ES13.3E3)') "Solndata vel  off proc after =", solndata(VELX_VAR:VELZ_VAR, neighInd(1),neighInd(2),neighInd(3))
!#endif
!        deallocate(xc)
!        deallocate(yc)
!        deallocate(zc)
      
!        call Grid_releaseBlkPtr(neighBlk, solndata)
        
!        call Eos_wrapped(MODE_DENS_EI, blkLimits, neighBlk)
        
!        messages = 0
        
!    end if
        
!  end do
    

!end if ! If Proc_ID = gr_meshMe

end subroutine inject_direct

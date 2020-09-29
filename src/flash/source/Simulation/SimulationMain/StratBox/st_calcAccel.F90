!!****if* source/physics/sourceTerms/Stir/StirMain/FromFile/st_calcAccel
!!
!! NAME
!!  st_calcAccel
!!
!! SYNOPSIS
!!  st_calcAccel(integer,INTENT (in)  :: blockid,
!!               integer,dimension(2,MDIM),INTENT(in)  :: blkLimitsGC,
!!               real(:),intent(IN) :: xcoord,
!!               real(:),intent(IN) :: ycoord,
!!               real(:),intent(IN) :: zcoord)
!!
!! DESCRIPTION
!!   Computes components of the driving force from Fourier modes
!!
!! ARGUMENTS
!!   blockID : ID of block in current processor
!!   blkLimitsGC : array holding the upper and lower index limits 
!!                 of an entire block (including GCs)
!!   xcoord : real array containg the coordinates along IAXIS
!!   ycoord : real array containg the coordinates along JAXIS
!!   zcoord : real array containg the coordinates along KAXIS
!!
!! AUTHOR
!!  Christoph Federrath, 2008
!!
!!***

subroutine st_calcAccel (blockID,blkLimitsGC,xCoord,yCoord,zCoord)

  use Stir_data, ONLY : st_maxmodes, st_mode, st_nmodes, st_aka, st_akb, st_solweightnorm, st_ampl
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr

  ! user modification: focus stirring to a thin layer in z.
  use Simulation_data, ONLY : sim_stirLayer, sim_stirH

  implicit none

#include "constants.h"
#include "Flash.h"

  ! Parameters

  integer,dimension(2,MDIM),INTENT(IN) :: blkLimitsGC
  integer,INTENT (IN)                  :: blockID

  real,dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),intent(IN) :: xCoord
  real,dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)),intent(IN) :: yCoord
  real,dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),intent(IN) :: zCoord

  ! pre-compute some mode info
  real,dimension(st_maxmodes,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) :: cosxi, sinxi
  real,dimension(st_maxmodes,blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) :: cosxj, sinxj
  real,dimension(st_maxmodes,blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: cosxk, sinxk

  real, pointer, dimension(:,:,:,:) :: solndata
  real                              :: realtrigterms, imtrigterms
  integer                           :: ib, ie
  integer                           :: i, j, k, m
  real, dimension(st_nmodes)        :: ampl

  ! user modification: focus stirring to a thin layer in z.
  real                              :: zz
  real                              :: layerWeight

  ! precompute some trig. (could be moved to Stir_init, but doesn't seem to save resources significantly)
  do m = 1, st_nmodes

     cosxk(m,:) = cos(st_mode(3,m)*zCoord(:))
     sinxk(m,:) = sin(st_mode(3,m)*zCoord(:))

     cosxj(m,:) = cos(st_mode(2,m)*yCoord(:))
     sinxj(m,:) = sin(st_mode(2,m)*yCoord(:))

     cosxi(m,:) = cos(st_mode(1,m)*xCoord(:))
     sinxi(m,:) = sin(st_mode(1,m)*xCoord(:))

     !! save some cpu time by precomputing one more multiplication
     ampl(m) = 2.0*st_ampl(m)

  enddo

  ib = blkLimitsGC(LOW,IAXIS)
  ie = blkLimitsGC(HIGH,IAXIS)

  call Grid_getBlkPtr(blockID,solnData)

  !! loop over all cells
  do k = blkLimitsGC (LOW, KAXIS), blkLimitsGC(HIGH,KAXIS)

     ! user modification: focus stirring to a thin layer in z.
     if (sim_stirLayer) then
        zz = zCoord(k)
        layerWeight = 1 + tanh((sim_stirH-abs(zz))/(0.5*sim_stirH))
     endif

     do j = blkLimitsGC (LOW, JAXIS), blkLimitsGC(HIGH,JAXIS)

#ifdef ACCX_VAR
        solnData(ACCX_VAR,ib:ie,j,k) = 0.0
#endif
#ifdef ACCY_VAR
        solnData(ACCY_VAR,ib:ie,j,k) = 0.0
#endif
#ifdef ACCZ_VAR
        solnData(ACCZ_VAR,ib:ie,j,k) = 0.0
#endif

        do i = blkLimitsGC (LOW, IAXIS), blkLimitsGC(HIGH,IAXIS)

           do m = 1, st_nmodes

                !  these are the real and imaginary parts, respectively, of
                !     e^{ i \vec{k} \cdot \vec{x} }
                !          = cos(kx*x + ky*y + kz*z) + i sin(kx*x + ky*y + kz*z)

                realtrigterms =    ( cosxi(m,i)*cosxj(m,j) - sinxi(m,i)*sinxj(m,j) ) * cosxk(m,k)  &
                                 - ( sinxi(m,i)*cosxj(m,j) + cosxi(m,i)*sinxj(m,j) ) * sinxk(m,k)

                imtrigterms   =    cosxi(m,i) * ( cosxj(m,j)*sinxk(m,k) + sinxj(m,j)*cosxk(m,k) ) &
                                 + sinxi(m,i) * ( cosxj(m,j)*cosxk(m,k) - sinxj(m,j)*sinxk(m,k) )
#ifdef ACCX_VAR
                solnData(ACCX_VAR,i,j,k) = solnData(ACCX_VAR,i,j,k) + ampl(m)*(st_aka(IAXIS,m)*realtrigterms-st_akb(IAXIS,m)*imtrigterms)
#endif
#ifdef ACCY_VAR
                solnData(ACCY_VAR,i,j,k) = solnData(ACCY_VAR,i,j,k) + ampl(m)*(st_aka(JAXIS,m)*realtrigterms-st_akb(JAXIS,m)*imtrigterms)
#endif
#ifdef ACCZ_VAR
                ! user modification: disable all vertical forcing
                !solnData(ACCZ_VAR,i,j,k) = solnData(ACCZ_VAR,i,j,k) + ampl(m)*(st_aka(KAXIS,m)*realtrigterms-st_akb(KAXIS,m)*imtrigterms)
#endif

           enddo  ! end loop over modes

#ifdef ACCX_VAR
           solnData(ACCX_VAR,i,j,k) = st_solweightnorm*solnData(ACCX_VAR,i,j,k)
           ! user modification: focus stirring to a thin layer in z.
           if (sim_stirLayer) then
              solnData(ACCX_VAR,i,j,k) = layerWeight*solnData(ACCX_VAR,i,j,k)
           endif
#endif
#ifdef ACCY_VAR
           solnData(ACCY_VAR,i,j,k) = st_solweightnorm*solnData(ACCY_VAR,i,j,k)
           ! user modification: focus stirring to a thin layer in z.
           if (sim_stirLayer) then
              solnData(ACCY_VAR,i,j,k) = layerWeight*solnData(ACCY_VAR,i,j,k)
           endif
#endif
#ifdef ACCZ_VAR
           ! user modification: disable all vertical forcing
           !solnData(ACCZ_VAR,i,j,k) = st_solweightnorm*solnData(ACCZ_VAR,i,j,k)
           ! user modification: focus stirring to a thin layer in z.
           !if (sim_stirLayer) then
           !   solnData(ACCZ_VAR,i,j,k) = layerWeight*solnData(ACCZ_VAR,i,j,k)
           !endif
#endif

        enddo  ! end loop over i
     enddo  ! end loop over j
  enddo  ! end loop over k

  call Grid_releaseBlkPtr(blockID,solnData)

  return

end subroutine st_calcAccel

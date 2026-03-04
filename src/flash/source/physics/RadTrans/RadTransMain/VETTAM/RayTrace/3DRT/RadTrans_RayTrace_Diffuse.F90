!===============================================================================
!
! Subroutine: RadTrans_raytrace_Diffuse
! Path: source/physics/RadTrans/RadTransMain/VET/RayTrace/sinkRT/RadTrans_raytrace_Diffuse.F90
!
!===============================================================================
!
! Description: Ray-Trace step for sink particles used with the VET cooling module. 
! Original Author: Manuel Jung (2018)
! Modified by Shyam Harimohan Menon to incorporate with VET module (2020-2021)
! Email: shyam.menon@anu.edu.au
!
!===============================================================================
!
! TO DO:
!
!===============================================================================
!
!===============================================================================

#include "Flash.h"
#include "constants.h"
subroutine RadTrans_raytrace_Diffuse()

  use RadTrans_RayTrace_3DRT
  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use calc_local_mod, ONLY : calc_local_block_contributions_3DRT
  use create_cut_block_mod, ONLY : create_cut_block_list_3DRT
  use Grid_data, ONLY: gr_oneBlock, gr_nBlockX, gr_nBlockY, gr_nBlockZ
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use rt_data_raytrace_3drt, ONLY : domainSizeX, domainSizeY, domainSizeZ, &
       nrOfAnglesPerGroup, nrOfAngles, dOmega
#ifdef XXED_VAR
  use RadTrans_data, ONLY: rt_etens
#endif
  use RadTrans_hybridCharModule

  implicit none
  !
  ! The number of coordinates to define one block.
  !
  integer, parameter  :: nrOfBlockCoords = SUM(2*KD)
  ! Some values needed for the solid angle integration of the
  ! radiation field. Phi and theta are the azimuthal and
  ! the polar angle respectively. They give the direction
  ! of the current ray, also called "characteristic".
  ! iAngle is the loop counter for the angle loop, which is
  ! embedded in the angle group loop. each angle group will
  ! contain the contribution from the angles in the index range
  ! from iAb to iAe. iAindex is used to store the face value
  ! contribution and regain them later.
  ! rt_nPhi and nTheta are the total numbers of azimuthal (phi)
  ! and polar angles (theta) respectively, which are used for the 
  ! solid angle integration of the radiation field. 
  !

  real                :: phi,theta
  integer             :: iAngle,iAb,iAe,iAindex
  real                :: dirFact(NDIM)

  ! Values for the use of solid angle groups. The number of
  ! angle groups defines for how many angles the face values
  ! are computed before they are communicated with an
  ! mpi_allgather command. Less communication is faster,
  ! but is also very memory consuming. One has to keep this
  ! in mind and choose the number of angle groups very carefully.
  !
  integer             :: iAngleGroup
  !
  !
  ! Method of Lambda-iteration:
  !     ALI = 0 | no acceleration.
  !     ALI = 1 | local Lambda-operator.
  !     ALI = 2 | tri-diagonal Lambda-operator (not-implemented).
  !
  ! The actual lambda operator, only local yet.
  !
  real                :: lambda_akku
  !
  ! The regular Factors are used to convert coordinates from
  ! physical values [cm] to cell units on the maximum level of
  ! refinement. The unit of regFact is then [cells cm^-1].
  ! Regular coordinates mean, that they are given in units
  ! of cell numbers on the maximum refinement level.
  !
  real,    save        :: regFact(NDIM)
  real,    save        :: PointSource(NDIM)

  ! The domain size in units of cells at the maximum
  ! refinement level.
  !
  integer, save        :: regTotI,regTotJ,regTotK

  !
  ! block indices
  !
  !    b: the index of blk in leafBlockList
  !  blk: the index of a block in the flash database
  !
  integer              :: b, blk
  !
  ! Some general flash loop iteration counter and parameters.
  !
  integer                         :: idx(NDIM), i, j, k, n, p, q
  real, dimension(MAXCELLS,NDIM) :: xr
  real :: taud, taud_local, taud_global
  real :: intensity, intensity_local, intensity_global

  !
  ! Physical constants from the database.
  !
  !
  ! All the rest
  !

  real :: dest(NDIM)
  !
  !-------------------------------------------------------------------------------
  !  
  !
  !======================================================================================
  !
  call Timers_start("rad_raytrace_diffuse")
  !
  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Starting Raytracing", "[3DRT] ")

  ! If we are doing iteration and we are at the first iteration step, 
  ! we have to contruct the lambda operator.
  !
  !if(iterationStep.eq.1 .and. rt_ALI.ne.0) then
  if(rt_ALI.ne.0) then
    do_lambda = .true.
  else
    do_lambda = .false.
  endif

  if(nrOfAngles .le. 0) return

  ! Initialize Angle Indices
  iAb=0
  iAe=0

  ! Get the current resolution on the highest refinement level
  ! and the conversion factors
  regTotI = gr_nBlockX * 2**(maxLevel-1) * NXB
  regTotJ = gr_nBlockY * 2**(maxLevel-1) * NYB
  regTotK = gr_nBlockZ * 2**(maxLevel-1) * NZB
  regFact(1) = real(regTotI) / domainSizeX
  regFact(2) = real(regTotJ) / domainSizeY
  regFact(3) = real(regTotK) / domainSizeZ
  !
  !================================================================================
  !
  ! All administrative stuff is done.
  ! We start with the loop over all angle groups.
  ! An angle group contains a loop over a specific range
  ! (nrOfAnglesPerAngleGroup) an integrate the specific intensities for these
  ! directions
  !
  !--------------------------------------------------------------------------------

  !Starting diffuse radiation contribution
  if(nrOfAngles.gt.0) then
    if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Background rays", "[3DRT] ")
    do iAngleGroup=1,rt_nrOfAngleGroups
      !
      if(dr_globalMe.eq.MASTER_PE) then
        write(io,*) ''
        write(io,*) 'Angle Group:',iAngleGroup
      endif
      !
#ifdef DEBUG_RT
      write(*,*) '--------------Angular Loop----------------------'
      write(*,*) 'rt_nPhi,rt_nTheta',rt_nPhi,rt_nTheta
      write(*,*) 'nrOfAngles',nrOfAngles
      write(*,*) 'dOmega',dOmega
#endif
      !
      ! Choose the angle range we are handling in this angle group
      !
      iAb = iAe+1 ! Start where we stopped 
      iAe = iAngleGroup*nrOfAnglesPerGroup
      !
#ifdef DEBUG_RT
      write(*,*) 'iAb,iAe',iAb,iAe
      write(*,*) '------------------------------------------------'
#endif

      ! Step 1 : Find the intensity (diffuse) and optical depth to faces of the block
      ! This would be stored, and used for the block-level solve later below
      call Timers_start("rad_calc_face_values")
      call rad_facevalues()
      call Timers_stop("rad_calc_face_values")

      if(dr_globalMe.eq.MASTER_PE) then
         write (io,*) 'Local Angle Loop finished'
      endif

      ! Step 2:  Gather and communicate the quantities at the faces to all processors
      call Timers_start("rad_MPI_ALLGATHER_faceValues")
      call gather_faces(iAb,iAe)
      call Timers_stop("rad_MPI_ALLGATHER_faceValues")


#ifdef DEBUG_RT
      write(*,*) 'starting Angular Loop for Global Contributions'
#endif

      ! Step 3: Calculate net contribution - local and global (i.e. from neighbouring faces)
      ! They are both done within the subroutine, including boundary conditions
      call Timers_start("rad_calcContrib")
      call rad_calcContrib()
      call Timers_stop("rad_calcContrib")

      if(dr_globalMe.eq.MASTER_PE) then
        write (io,*) 'Global Angle Loop finished'
        write (io,*)
      endif
      !
   enddo ! iAngleGroup
 end if

  ! Calculate Variable Eddington Tensor
#ifdef XXED_VAR
  if(rt_etens .eq. "vet") then
    if(nrOfAngles.gt.0) then
      do b=1, nrOfLeafBlocks
        blk=leafList(b)
        call Grid_getBlkPtr(blk, solnData) 
        do k=ib(3)+K3D,ie(3) ! NOTE: kb+1 for 3D
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)

              if(solnData(MEAN_VAR,i,j,k).GT.0.) then 
                solnData(XXED_VAR,i,j,k) = solnData(RPXX_VAR,i,j,k)/solnData(MEAN_VAR,i,j,k)
#if NDIM>1
                solnData(XYED_VAR,i,j,k) = solnData(RPXY_VAR,i,j,k)/solnData(MEAN_VAR,i,j,k)
                solnData(YYED_VAR,i,j,k) = solnData(RPYY_VAR,i,j,k)/solnData(MEAN_VAR,i,j,k)
#if NDIM>2
                solnData(XZED_VAR,i,j,k) = solnData(RPXZ_VAR,i,j,k)/solnData(MEAN_VAR,i,j,k)
                solnData(YZED_VAR,i,j,k) = solnData(RPYZ_VAR,i,j,k)/solnData(MEAN_VAR,i,j,k)
                solnData(ZZED_VAR,i,j,k) = solnData(RPZZ_VAR,i,j,k)/solnData(MEAN_VAR,i,j,k)
#endif
#endif
              !If MEAN intensity is l.t. 0 something is wrong. Exit...
              else if(solnData(MEAN_VAR,i,j,k).LT.0.) then
                call Driver_abortFlash("MEAN_VAR value found to be less than zero. Exiting...")

              !MEAN_VAR is zero here, so will be the 2nd moment. Set f=f_Edd in this case
              else

                solnData(XXED_VAR,i,j,k) = 1./3.
#if NDIM>1
                solnData(XYED_VAR,i,j,k) = 0.0
                solnData(YYED_VAR,i,j,k) = 1./3.
#if NDIM>2
                solnData(XZED_VAR,i,j,k) = 0.0
                solnData(YZED_VAR,i,j,k) = 0.0
                solnData(ZZED_VAR,i,j,k) = 1./3.
#endif
#endif

              endif
            end do 
          end do
        end do
        call Grid_releaseBlkPtr(blk, solnData)
      end do
    end if
  end if

#endif

  if(dr_globalMe.eq.MASTER_PE) then
     write(io,*)
     write(io,*) 'raytrace_3DRT done'
     write(io,*)
     write(io,*) '======End of raytrace_3DRT======'
     write(io,*)
     write(io,*)
  endif
  !
  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Finished Raytracing", "[3DRT] ")
  !
  call Timers_stop("rad_raytrace_diffuse")
  
  CONTAINS

!
! First: Calculate all face values for the intensity (meaf) and
! optical depth (tauf) for all angles in this angle group
  SUBROUTINE rad_facevalues()
    implicit none
    integer :: iter
    real :: dtaucell, dlcell
       !
       !    Store (i.e. pack) the local face values and block coordinates.
       !
#ifdef DEBUG_RT
       write(*,*) 'storing face values for communication'
#endif

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)

      do n=IAXIS,KAXIS
       call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do
      do iAngle=iAb,iAe
       iAindex = (iAngle-iAb)+1
       call getAngles(iAngle,rand_angles(1),rand_angles(2),rand_angles(3),theta,phi,dirFact)
       do n=1,NDIM
         p = MOD(n,3) + 1
         q = MOD(n+1,3) + 1

         IF(dirFact(n).lt.0.) THEN
           i = ib(n)
         ELSE
           i = ie(n)
         END IF
         do j = ib(q),ie(q)
           do k = ib(p),ie(p)
             idx(n) = i
             idx(q) = j
             idx(p) = k
             ! Here is space for optimization
             dest = (/ xr(idx(1),1), xr(idx(2),2), xr(idx(3),3) /)
             call calc_local_block_contributions_3DRT(    &
                  regFact, blk, &
                  PointSource, dest, dirFact,   &
                  maxLevel,         &
                  taud, intensity, lambda_akku, &
                  n, &
                  .false.,dtaucell,dlcell)
             faceValueAll(k,j,n,1,b,iAindex,dr_globalMe) = taud
             faceValueAll(k,j,n,2,b,iAindex,dr_globalMe) = intensity
           end do ! k
         end do ! j
       end do ! n
      end do ! iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    end do ! b
#ifdef DEBUG_RT
       write(*,*) 'storing face values for communication finished'
#endif
  end SUBROUTINE rad_facevalues

  SUBROUTINE rad_calcContrib
    use SemenovOpacities, ONLY: getOpacity_planck
    implicit none
    real :: dtaucell, dlcell
       !
#ifdef DEBUG_RT
       write(*,*) 'calculating cut block contributions'
       write(*,*) 'iAngle:',iAngle
#endif
       !
   do b=1,nrOfLeafBlocks
     blk = leafList(b)
     call Grid_getBlkPtr(blk, solnData)

     do iAngle=iAb,iAe

       iAindex = (iAngle-iAb)+1
       call getAngles(iAngle,rand_angles(1),rand_angles(2),rand_angles(3),&
        theta,phi,dirFact)

       do k=ib(3)+K3D,ie(3) ! NOTE: kb+1 for 3D
         do j=ib(2)+K2D,ie(2)
           do i=ib(1)+1,ie(1)
             dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
             dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
             dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)
             !
             ! Second: Calculate all local contributions to the intensity (mean) and
             ! optical depth (taud) for all angles in this angle group.
             !

             call calc_local_block_contributions_3DRT(    &
                  regFact, blk, &
                  PointSource, dest, dirFact, &
                  maxLevel,                               &
                  taud_local, intensity_local, lambda_akku, &
                  0, .false.,dtaucell,dlcell)


             !
             ! Raytrace through the block structure and add
             ! contributions from the corner (Lars's face)
             ! values
             ! Store the result in taud_global, intensity_global.

             ! irradiation is added below and extinct with the final optical depth to the destination
             call rt_irradiation(dirFact, intensity_global)

             call create_cut_block_list_3DRT(                  &
                  dr_globalMe, blk,                            &
                  taud_global, intensity_global,               &
                  PointSource, dest, &
                  regFact, dirFact, &
                  iAindex, &
                  maxNrOfLeafBlocks, &
                  dr_globalNumProcs, maxLevel, &
                  .false.)

             ! Now add up global and local optical depth
             taud = taud_global + taud_local

             ! Now add up global and local intensities.
             ! intensity now contains the specific intensity for the current
             ! direction.
             !-------------------------------------
             intensity = intensity_local + intensity_global * exp(-taud_local)

             ! Add specific intensity (intensity) to the mean intensity integral
             solnData(MEAN_VAR,i,j,k) = solnData(MEAN_VAR,i,j,k) &
               + MAX(1.0/(4.0*PI) * intensity * dOmega,0.)

             ! Add to 2nd moments for VET calculation
             solnData(RPXX_VAR,i,j,k) = solnData(RPXX_VAR,i,j,k) &
             + 1.0/(4.0*PI)*intensity*dirFact(IAXIS)*dirFact(IAXIS)*dOmega
#if NDIM>1
             solnData(RPXY_VAR,i,j,k) = solnData(RPXY_VAR,i,j,k) &
             + 1.0/(4.0*PI)*intensity*dirFact(IAXIS)*dirFact(JAXIS)*dOmega
             solnData(RPYY_VAR,i,j,k) = solnData(RPYY_VAR,i,j,k) &
             + 1.0/(4.0*PI)*intensity*dirFact(JAXIS)*dirFact(JAXIS)*dOmega  
#if NDIM>2
             solnData(RPXZ_VAR,i,j,k) = solnData(RPXZ_VAR,i,j,k) &
             + 1.0/(4.0*PI)*intensity*dirFact(IAXIS)*dirFact(KAXIS)*dOmega
             solnData(RPYZ_VAR,i,j,k) = solnData(RPYZ_VAR,i,j,k) &
             + 1.0/(4.0*PI)*intensity*dirFact(JAXIS)*dirFact(KAXIS)*dOmega
             solnData(RPZZ_VAR,i,j,k) = solnData(RPZZ_VAR,i,j,k) &
             + 1.0/(4.0*PI)*intensity*dirFact(KAXIS)*dirFact(KAXIS)*dOmega
#endif
#endif

#ifdef LAMB_VAR
             if(do_lambda) then
               solnData(LAMB_VAR,i,j,k) = solnData(LAMB_VAR,i,j,k) &
                 + lambda_akku * 1.d0/(4.d0*PI) * dOmega
             endif
#endif
            enddo ! i
          enddo ! j
        enddo ! k
      enddo ! iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! blocks
    
#ifdef DEBUG_RT
    write(*,*) 'calculating cut block contributions finished'
#endif
  end SUBROUTINE rad_calcContrib

  SUBROUTINE gather_faces(iAb,iAe)
    implicit none
    integer, intent(in) :: iAb,iAe
    integer :: angles
    integer, dimension(dr_globalNumProcs) :: displs,recvcounts
    !
    !  Gather all face values from all processors and store
    !  the results in faceValueAll.
    !
#ifdef DEBUG_RT
    write(*,*) 'gathering face-values...'
#endif
    angles = iAe-iAb+1
    IF(angles.eq.nrOfAnglesPerGroup) THEN
      call MPI_ALLGATHER( &
           MPI_IN_PLACE, &
           0, &
           MPI_DATATYPE_NULL, &
           faceValueAll,&
           (nb(1)+1)*(nb(2)+1)*3*2*maxNrOfLeafBlocks*nrOfAnglesPerGroup, &
           FLASH_REAL,&
           MPI_COMM_WORLD,ierr)
    ELSE IF(iAe-iAb+1.ge.1) THEN
      do i=1,dr_globalNumProcs
        displs(i) = (nb(1)+1)*(nb(2)+1)*3*2*maxNrOfLeafBlocks*nrOfAnglesPerGroup*(i-1)
      end do
      recvcounts(:) = (nb(1)+1)*(nb(2)+1)*3*2*maxNrOfLeafBlocks*angles

      call MPI_ALLGATHERV( &
           MPI_IN_PLACE, &
           0, &
           MPI_DATATYPE_NULL, &
           faceValueAll, &
           recvcounts, &
           displs, &
           FLASH_REAL, &
           MPI_COMM_WORLD,ierr)
    END IF
#ifdef DEBUG_RT
    write(*,*) 'gathering face-values... done'
#endif
    !
  END SUBROUTINE gather_faces


end subroutine RadTrans_Raytrace_Diffuse
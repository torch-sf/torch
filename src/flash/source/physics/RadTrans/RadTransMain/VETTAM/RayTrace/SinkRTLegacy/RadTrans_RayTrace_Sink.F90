!===============================================================================
!
! Subroutine: RadTrans_raytrace_Sink
! Path: source/physics/RadTrans/RadTransMain/VET/RayTrace/sinkRT/RadTrans_raytrace_Sink.F90
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
!

#undef RT_CORRECT
#undef DEBUG_RT

#include "Flash.h"
#include "constants.h"
subroutine RadTrans_raytrace_Sink()
  use RadTrans_RayTrace_3DRT
  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use calc_local_sink_mod, ONLY: calc_local_block_contributions_sink
  use create_cut_block_mod, ONLY : create_cut_block_list_3DRT
  use Grid_data, ONLY: gr_oneBlock, gr_nBlockX, gr_nBlockY, gr_nBlockZ
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use rt_data_raytrace_3drt, ONLY : domainSizeX, domainSizeY, domainSizeZ, &
       nrOfAnglesPerGroup, nrOfAngles
  use Particles_sinkData, ONLY : particles_global, localnp, localnpf
  use pt_sinkSort,      ONLY: NewQsort_IN
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal
  use RadTrans_data, ONLY: rt_speedlt
  use RadTrans_hybridCharModule

  implicit none

  integer             :: iAngle,iAb,iAe,iAindex,iAnglemax
  real                :: dirFact(NDIM)
  integer             :: iAngleGroup
  ! The regular Factors are used to convert coordinates from
  ! physical values [cm] to cell units on the maximum level of
  ! refinement. The unit of regFact is then [cells cm^-1].
  ! Regular coordinates mean, that they are given in units
  ! of cell numbers on the maximum refinement level.
  !
  real,    save       :: regFact(NDIM)
  real,    save       :: PointSource(NDIM)

  ! The domain size in units of cells at the maximum
  ! refinement level.
  !
  integer, save       :: regTotI,regTotJ,regTotK

  !
  ! block indices
  !
  !  b: the index of blk in leafBlockList
  !  blk: the index of a block in the flash database
  !
  integer             :: b, blk

  ! The actual lambda operator, only local yet.
  !
  real                :: lambda_akku

  !
  ! Some general flash loop iteration counter and parameters.
  !
  integer             :: idx(NDIM), i, j, k, n, p, q
  real                :: distSqr_sink
  real, dimension(MAXCELLS,NDIM) :: xr
  real :: taud, taud_local, taud_global
  real :: intensity, intensity_local, intensity_global
  real :: dest(NDIM)

  integer, dimension(:), allocatable :: id_sorted, QSindex
  ! Particle stuff
  integer              :: nAllParticles



  !-------------------------------------------------------------------------------
  !  
  !
  !======================================================================================
  !

  ! gather all particle positions and luminosities
  !
#ifdef SINK_PART_TYPE
  call pt_sinkGatherGlobal()
  nAllParticles = localnpf
#else
  nAllParticles = 0
#endif

  !Return if no sink particles present
  if(nAllParticles .eq. 0) return

  call Timers_start("rad_raytrace_sink")

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Starting Raytracing for sinks", "[RT_SINK] ")

  call rad_resetsinkContainers()

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

  ! particles_global is not in the same order on all procs. Order it!
  allocate(id_sorted(localnpf),QSindex(localnpf))
  do iAngle=1,localnpf
    id_sorted(iAngle) = int(particles_global(TAG_PART_PROP,iAngle))
  end do
  call NewQsort_IN(id_sorted, QSindex)

  !Make sure the communication is limited to the memory-related constraint
  iAnglemax = CEILING(nAllParticles/REAL(nrOfAnglesPerGroup))
  do iAngleGroup=1,iAnglemax
    iAb = (iAngleGroup-1)*nrOfAnglesPerGroup + 1
    iAe = MIN(iAngleGroup*nrOfAnglesPerGroup,nAllParticles)

    call Timers_start("rad_calc_face_values")
    call rad_facevalues()
    call Timers_stop("rad_calc_face_values")

    call Timers_start("rad_MPI_ALLGATHER_faceValues")
    call gather_faces(iAb,iAe)
    call Timers_stop("rad_MPI_ALLGATHER_faceValues")

    call Timers_start("rad_calcContrib")
    call rad_calcContrib()
    call Timers_stop("rad_calcContrib")
  end do

  deallocate(id_sorted,QSindex)

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Raytracing for sinks completed", "[RT_SINK] ")

  call Timers_stop("rad_raytrace_sink")

  CONTAINS

  subroutine rad_facevalues()
    real :: dtaucell, dlcell

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      do n=IAXIS,KAXIS
        call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do

      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        
        !Position of point source
        PointSource = particles_global(POSX_PART_PROP:POSZ_PART_PROP,QSindex(iAngle))
        
        !Set opacity for this sink particle
        call Timers_start("setSinkOpacity")
        call SetSinkOpacity(blk,QSindex(iAngle))
        call Timers_stop("setSinkOpacity")

        ! Compute optical depths to faces
        ! The optical depth should be zero, if either:
        ! - the destination is equal the source
        ! - the face can not be reached in that direction
        ! (in all blocks but the on including the point source,
        !  only three faces have valid taud values.)
        do n=1,NDIM
          !Direction of face required. For x -y,z, for y: z,x, for z: x,y
          p = MOD(n,3) + 1
          q = MOD(n+1,3) + 1

          i=ib(n)
          do j = ib(q),ie(q)
            do k = ib(p),ie(p)
              idx(n) = i
              idx(q) = j
              idx(p) = k
              ! Here is space for optimization
              dest = (/ xr(idx(1),1), xr(idx(2),2), xr(idx(3),3) /)
              dirFact = dest-PointSource
              
              if(dirFact(n) .ge. 0) then
                taud = 0.0
              else
                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, taud, &
                        n,dtaucell,dlcell,.false.)
              end if
              faceValueAll(k,j,n,1,b,iAindex,dr_globalMe) = taud
            end do ! k
          end do !j


          i=ie(n)
          do j = ib(q),ie(q)
            do k = ib(p),ie(p)
              idx(n) = i
              idx(q) = j
              idx(p) = k
              ! Here is space for optimization
              dest = (/ xr(idx(1),1), xr(idx(2),2), xr(idx(3),3) /)
              dirFact = dest-PointSource
              if(dirFact(n) .lt. 0) then
                taud = 0.0
              else
                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, taud, &
                        n,dtaucell,dlcell,.false.)
              end if
              faceValueAll(k,j,n,2,b,iAindex,dr_globalMe) = taud
            end do ! k
          end do !j

        end do !n
      end do !iAngle
    enddo ! b

  END SUBROUTINE rad_facevalues

  SUBROUTINE rad_calcContrib
    real :: luminosity, localtaud, rho, del(NDIM), radius_sink, flux
    real :: opac, l_segment
    real :: phi,theta
    real :: dtaucell, dlcell
    logical :: printval
    REAL,DIMENSION(NDIM) :: target_point

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      call Grid_getDeltas(blk, del)
      do n=IAXIS,KAXIS
        call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do

      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        
        !Position of point source
        PointSource = particles_global(POSX_PART_PROP:POSZ_PART_PROP,QSindex(iAngle))
        
        !Set opacity for this sink particle
        call Timers_start("setSinkOpacity")
        call SetSinkOpacity(blk,QSindex(iAngle))
        call Timers_stop("setSinkOpacity")

        do k=ib(3)+K3D,ie(3)
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)
              dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
              dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
              dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)

              dirFact = dest-PointSource
              dirFact = dirFact/SQRT(SUM(dirFact**2))

              !DEBUG Stuff
#ifdef DEBUG_RT
              !Choose target point for which to print quantities
              target_point = (/del(1)/2.,-del(2)/2.,-del(3)/2./)
              if((target_point(1) .gt. dest(1) - del(1)/2. .and. target_point(1) .lt. dest(1) + del(1)/2.) &
                  .and. (target_point(2) .gt. dest(2) - del(2)/2. .and. target_point(2) .lt. dest(2) + del(2)/2.) & 
                  .and. (target_point(3) .gt. dest(3) - del(3)/2. .and. target_point(3) .lt. dest(3) + del(3)/2.)) then
                printval = .true.
                print *, 'Debugging and printing info for point:',dest
              else
                printval = .false.
              endif
#else
              printval = .false.
#endif


              call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, taud_local, &
                        0,dtaucell,dlcell,printval)

              call create_cut_block_list_3DRT(&
                        dr_globalMe, blk,taud_global, intensity_global, PointSource, dest, &
                        regFact, dirFact, iAindex, maxNrOfLeafBlocks, dr_globalNumProcs, maxLevel, .true.)

              opac = solnData(OPAC_VAR,i,j,k)

              ! Now add up global and local optical depth
              taud = taud_global + taud_local

              ! Now add up global and local intensities.
              ! intensity now contains the specific intensity for the current
              ! direction.
              !-------------------------------------
              ! Take the maximum to limit the distance to bigger distances then 1./2 of a cell
              ! This is approximately the average, if we calculate the distance between any point
              ! in the cell and the star particle.

              ! SHM: This doesn't really come into play if the zero opacity condition within the accretion radius of the sink is used
              distSqr_sink = MAX(SUM((dest-PointSource)**2),(0.5/regFact(1))**2)

              !Prevent any radiation within the radius of the star
              radius_sink = particles_global(STELLAR_RADIUS_PART_PROP,QSindex(iAngle))
              if (distSqr_sink .ge. radius_sink**2) then

               intensity = particles_global(LUMINOSITY_PART_PROP,QSindex(iAngle)) / &
                (4.0*PI*distSqr_sink) * exp(-taud)

              else 

               intensity = 0.0

              end if

              flux = intensity
              !Use returned values from calc_local for these
              localtaud = dtaucell
              !Segment length halved to keep consistent with earlier versions of code
              l_segment = dlcell/2.
              IF(localtaud.lt.1.e-4) THEN
                flux = flux * opac
              ELSE
                flux = flux * (1.-exp(-localtaud))/(2.*l_segment)
              END IF

              !Radiation heating term due to stellar source
              solnData(FSPT_VAR,i,j,k) = solnData(FSPT_VAR,i,j,k) + flux
#ifdef DEBUG_RT
              if(printval) print * , 'taud_local, taud_global, localtaud, l_segment, fspt',taud_local,&
                taud_global, localtaud, l_segment, solnData(FSPT_VAR,i,j,k)
#endif

              !Radiation pressure source terms due to stellar source
              solnData(STMX_VAR,i,j,k) = solnData(STMX_VAR,i,j,k) + solnData(FSPT_VAR,i,j,k) /&
                                                rt_speedlt * dirFact(IAXIS)
#if NDIM>1
              solnData(STMY_VAR,i,j,k) = solnData(STMY_VAR,i,j,k) + solnData(FSPT_VAR,i,j,k) /&
                                                rt_speedlt * dirFact(JAXIS)
#if NDIM>2
              solnData(STMZ_VAR,i,j,k) = solnData(STMZ_VAR,i,j,k) + solnData(FSPT_VAR,i,j,k) /&
                                                rt_speedlt * dirFact(KAXIS)
#endif
#endif

              solnData(TAUD_VAR,i,j,k) = taud
              solnData(TAUF_VAR,i,j,k) = taud_global
              solnData(TAUZ_VAR,i,j,k) = taud_local

            enddo ! i
          enddo !j
        enddo !k
      enddo !iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! blocks

  END SUBROUTINE rad_calcContrib

  subroutine gather_faces(iAb,iAe)
    integer, intent(in) :: iAb,iAe
    integer :: angles
    integer, dimension(dr_globalNumProcs) :: displs,recvcounts
    !
    !  Gather all face values from all processors and store
    !  the results in faceValueAll.
    !
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
    !
  END SUBROUTINE gather_faces

  SUBROUTINE SetSinkOpacity(blk,pno)
    use Opacity_sink_interface, ONLY: getOpacity_draine
    use Particles_sinkData, ONLY: accretion_radius
    use SemenovOpacities, ONLY: getOpacity_star
    use RadTrans_data, ONLY: rt_stellarop_type,rt_stellar_opacity,rt_boltz
    IMPLICIT NONE
    integer, intent(in) :: blk, pno
    real :: rho, temp, kappa, Tstar, lum_star, rad_star, xdist, ydist, zdist, distsqr_point, dist_cell
    integer :: i,j,k
    real, dimension(NDIM) :: position, sink_pos

    if(rt_stellarop_type .eq. "Draine") then 
      !Get color temperature of star from luminosity and radius
      lum_star = particles_global(LUMINOSITY_PART_PROP,pno)
      rad_star = particles_global(STELLAR_RADIUS_PART_PROP,pno)
      Tstar = (lum_star/(4.*PI*rad_star**2*rt_boltz))**0.25
      call getOpacity_draine(Tstar,kappa)
    else if(rt_stellarop_type .eq. "Fixed") then
      kappa = rt_stellar_opacity
    else if(rt_stellarop_type .eq. "Semenov") then 
      continue
    else
      call Driver_abortFlash("RadTrans_raytrace_3DRT.F90: Stellar opacity type not recognised.")
    endif

    do k=GRID_KLO_GC,GRID_KHI_GC
      do j=GRID_JLO_GC,GRID_JHI_GC
        do i=GRID_ILO_GC,GRID_IHI_GC
          !Get distance from sink centre
          !Sink position
          sink_pos = particles_global(POSX_PART_PROP:POSZ_PART_PROP,pno)
          !Cell position
          position(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
          position(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
          position(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)
          !Distance
          xdist = position(1) - sink_pos(1)
          ydist = position(2) - sink_pos(2)
          zdist = position(3) - sink_pos(3)
          distsqr_point = xdist**2 + ydist**2 + zdist**2
          dist_cell = SQRT(distsqr_point)

          rho  = solnData(DENS_VAR,i,j,k)
          
          if(rt_stellarop_type .eq. "Semenov") then
#ifdef TEMP_VAR
            temp = solnData(TEMP_VAR,i,j,k)
#endif
            call getOpacity_star(pno,temp,rho,kappa)
          endif

          !Set cells within the accretion radius to be transparent
          if(dist_cell .gt. accretion_radius) then
            solnData(OPAC_VAR,i,j,k) = kappa * rho
          else
            solnData(OPAC_VAR,i,j,k) = 0.0
          endif

        end do 
      end do
    end do

END SUBROUTINE SetSinkOpacity

END SUBROUTINE RadTrans_raytrace_Sink





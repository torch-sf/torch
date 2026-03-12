!===============================================================================
!
! Subroutine: RayTrace_LW
! Path: source/physics/RadTrans/RadTransMain/VETTAM/LWBand/RayTrace_LW.F90
!
!===============================================================================
!
! Description: Ray-Trace step for computing the dissociation of H2 via LW photons. 
! Author: Shyam Menon (2023)
! Email: shyam.menon@rutgers.edu
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
subroutine RayTrace_LW()
  use RadTrans_RayTrace_3DRT
  use rt_lwdata
  use rt_lwmodule, ONLY: get_fshield_H, get_fshield_H2
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

  !
  ! Some general flash loop iteration counter and parameters.
  !
  integer             :: idx(NDIM), i, j, k, n, p, q
  real                :: distSqr_sink
  real, dimension(MAXCELLS,NDIM) :: xr
  real :: ColDens, ColDens_local, ColDens_global
  real :: intensity_global
  real :: dest(NDIM)

  integer, dimension(:), allocatable :: id_sorted, QSindex
  ! Particle stuff
  integer              :: nAllParticles
  !Flag for which species the shielding is computed; 0: H2, 1: H
  integer              :: fsh_spec



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

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Starting Raytracing for sink Lyman Werner contribution", "[RT_SINK] ")

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

  fsh_spec = 0 !Indicate that shielding function for H2 is to be done
  call rad_resetLWContainers() 
  !Set the variable OPAC_VAR to the volume number density
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call SetOPAC_LW(blk)
  end do
    

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

  !Compute the LW-flux weighted average shielding factor for H2
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do

  if(useHIshield) then
    !Now do the same steps for H
    fsh_spec = 1 !Indicate that shielding function for H is to be done
    call rad_resetLWContainers() ! This needs to be done to reset the LWFL_VAR container and the FSH
    !Set the variable OPAC_VAR to the volume number density
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call SetOPAC_LW(blk)
    end do
    
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

    !Compute the LW-flux weighted average shielding factor for H
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call average_fshield(blk)
    end do
  endif


  fsh_spec = 2 !This is done to ensure only LWFL_VAR is reset in next step
  call rad_resetLWContainers()
  !Now combine above shielding factors to get the LW flux
  call calculate_LWFlux()

  deallocate(id_sorted,QSindex)
    

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Raytracing for Lyman Werner Radiation completed", "[RT_SINK] ")

  call Timers_stop("rad_raytrace_sink")

  CONTAINS

  subroutine rad_facevalues()
    real :: dColDens, dlcell

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

        ! Compute optical depths to faces
        ! The optical depth should be zero, if either:
        ! - the destination is equal the source
        ! - the face can not be reached in that direction
        ! (in all blocks but the on including the point source,
        !  only three faces have valid ColDens values.)
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
                ColDens = 0.0
              else
                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, ColDens, &
                        n,dColDens,dlcell,.false.)
              end if
              faceValueAll(k,j,n,1,b,iAindex,dr_globalMe) = ColDens
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
                ColDens = 0.0
              else
                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, ColDens, &
                        n,dColDens,dlcell,.false.)
              end if
              faceValueAll(k,j,n,2,b,iAindex,dr_globalMe) = ColDens
            end do ! k
          end do !j

        end do !n
      end do !iAngle
    enddo ! b

  END SUBROUTINE rad_facevalues

  SUBROUTINE calculate_LWFlux
    real :: luminosity, del(NDIM), radius_sink, flux, fshield_eff

    !Loop over all blocks and sinks
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

        do k=ib(3)+K3D,ie(3)
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)
                dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
                dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
                dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)

                dirFact = dest-PointSource
                dirFact = dirFact/SQRT(SUM(dirFact**2))

                distSqr_sink = MAX(SUM((dest-PointSource)**2),(0.5/regFact(1))**2)

                !Prevent any radiation within the radius of the star
                radius_sink = particles_global(STELLAR_RADIUS_PART_PROP,QSindex(iAngle))
                if (distSqr_sink .ge. radius_sink**2) then
#ifdef LUMLW_PART_PROP
                    luminosity = particles_global(LUMLW_PART_PROP,QSindex(iAngle))
#elif NPEP_PART_PROP
                    luminosity = particles_global(NPEP_PART_PROP,QSindex(iAngle))
#else
                    luminosity = particles_global(LUMINOSITY_PART_PROP,QSindex(iAngle))
#endif
                
                    !Unattenuated radiation flux
                    flux = luminosity / (4.0*PI*distSqr_sink)

                else 
                    flux = 0.0

                end if

                fshield_eff = solnData(FSHM_VAR,i,j,k)
                !Include shielding due to atomic H if useHIShield included
                if(useHIshield) fshield_eff = fshield_eff * solnData(FSHA_VAR,i,j,k) 

                solnData(LWFL_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) + flux * fshield_eff
              enddo ! i
          enddo !j
        enddo !k
      enddo !iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! blocks
  END SUBROUTINE calculate_LWFlux



  SUBROUTINE rad_calcContrib
    real :: luminosity, localColDens, rho, del(NDIM), radius_sink, flux, fshield
    real :: dColDens, dlcell
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

        do k=ib(3)+K3D,ie(3)
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)
                dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
                dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
                dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)

                dirFact = dest-PointSource
                dirFact = dirFact/SQRT(SUM(dirFact**2))

                printval = .false.

                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, ColDens_local, &
                        0,dColDens,dlcell,printval)

                call create_cut_block_list_3DRT(&
                        dr_globalMe, blk,ColDens_global, intensity_global, PointSource, dest, &
                        regFact, dirFact, iAindex, maxNrOfLeafBlocks, dr_globalNumProcs, maxLevel, .true.)

                ! Now add up global and local column densities
                ColDens = ColDens_global + ColDens_local

                ! Now add up global and local column densities.
                !-------------------------------------
                ! Take the maximum to limit the distance to bigger distances then 1./2 of a cell
                ! This is approximately the average, if we calculate the distance between any point
                ! in the cell and the star particle.

                distSqr_sink = MAX(SUM((dest-PointSource)**2),(0.5/regFact(1))**2)

                !Prevent any radiation within the radius of the star
                radius_sink = particles_global(STELLAR_RADIUS_PART_PROP,QSindex(iAngle))
                if (distSqr_sink .ge. radius_sink**2) then
#ifdef LUMLW_PART_PROP
                    luminosity = particles_global(LUMLW_PART_PROP,QSindex(iAngle))
#elif NPEP_PART_PROP
                    luminosity = particles_global(NPEP_PART_PROP,QSindex(iAngle))
#else
                    luminosity = particles_global(LUMINOSITY_PART_PROP,QSindex(iAngle))
#endif
                
                    !Unattenuated radiation flux
                    flux = luminosity / (4.0*PI*distSqr_sink)

                else 
                    flux = 0.0

                end if

                !Get the shielding factor for the obtained column density
                if(fsh_spec == 0) then
                  fshield = get_fshield_H2(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(FSHM_VAR,i,j,k) = solnData(FSHM_VAR,i,j,k) + flux
                else if(fsh_spec == 1) then
                  fshield = get_fshield_H(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(FSHA_VAR,i,j,k) = solnData(FSHA_VAR,i,j,k) + flux
                endif

                !Flux-weighted shielding: this is not the final value of the LW flux used
                !This counter is used to compute the average flux-weighted fshield over all species; will be replaced later
                solnData(LWFL_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) + flux * fshield
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

  SUBROUTINE SetOPAC_LW(blk)
    use rt_ionisedata, ONLY: h2A, hA
    IMPLICIT NONE
    integer, intent(in) :: blk
    real :: nH2, nH
    integer :: i,j,k
    real, dimension(NDIM) :: position, sink_pos

    call Grid_getBlkPtr(blk, solnData)
    do k=GRID_KLO_GC,GRID_KHI_GC
      do j=GRID_JLO_GC,GRID_JHI_GC
        do i=GRID_ILO_GC,GRID_IHI_GC
            
            if(fsh_spec == 0) then
              !Obtain the number density of nH2
#ifdef H2_SPEC
              nH2 = solnData(DENS_VAR,i,j,k) * solnData(H2_SPEC,i,j,k)/h2A
#else
              call Driver_abortFlash("The Multispecies contained H2_SPEC has to be present for LW dissociation.")
#endif
              !Now set the container OPAC_VAR which is the quantity integrated in the ray-tracer to the number density;
              !this ensures the column density is returned
              solnData(OPAC_VAR,i,j,k) = nH2
            else if(fsh_spec == 1) then
              !Obtain the number density of nH
#ifdef IHA_SPEC
              nH = solnData(DENS_VAR,i,j,k) * solnData(IHA_SPEC,i,j,k)/hA
#else
              call Driver_abortFlash("The Multispecies contained IHA_SPEC has to be present for LW dissociation.")
#endif
              !Now set the container OPAC_VAR which is the quantity integrated in the ray-tracer to the number density;
              !this ensures the column density is returned
              solnData(OPAC_VAR,i,j,k) = nH
            endif

        end do 
      end do
    end do
    call Grid_releaseBlkPtr(blk, solnData)

  END SUBROUTINE SetOPAC_LW

  SUBROUTINE average_fshield(blk)
    IMPLICIT NONE
    integer, intent(in) :: blk
    integer :: i,j,k
    real :: totalflux, fshield

    call Grid_getBlkPtr(blk, solnData)
    do k=GRID_KLO_GC,GRID_KHI_GC
        do j=GRID_JLO_GC,GRID_JHI_GC
            do i=GRID_ILO_GC,GRID_IHI_GC
                !The total cumulative flux from sinks is stored in FSHM_VAR or FSHA_VAR (depending on which is being updated now)
                if(fsh_spec == 0) then
                  totalflux = solnData(FSHM_VAR,i,j,k)
                else if(fsh_spec == 1) then
                  totalflux = solnData(FSHA_VAR,i,j,k)
                endif
                !Now replace this var with the appropriate data, i.e. the avg shielding factor in FSHL_VAR
                if(totalflux .gt. 0) then
                    fshield = solnData(LWFL_VAR,i,j,k)/totalflux
                else
                    fshield = 0.0
                endif
                
                !Now store the shielding quantities in their containers
                if(fsh_spec == 0) then
                  solnData(FSHM_VAR,i,j,k) = fshield
                else if(fsh_spec == 1) then
                  solnData(FSHA_VAR,i,j,k) = fshield
                endif
                    
            end do
        end do
    end do
    call Grid_releaseBlkPtr(blk, solnData)

  END SUBROUTINE average_fshield

  SUBROUTINE rad_resetLWContainers
    implicit none
    INTEGER :: b, blk

    call Timers_start("rad_resetData")

    do b=1,nrOfLeafBlocks
      !
      blk = leafList(b)

      call Grid_getBlkPtr(blk, solnData)
      !Molecular
      if(fsh_spec == 0) solnData(FSHM_VAR,:,:,:) = 0.0
      if(fsh_spec == 1) solnData(FSHA_VAR,:,:,:) = 0.0
      solnData(LWFL_VAR,:,:,:) = 0.0
      call Grid_releaseBlkPtr(blk, solnData)

    enddo ! b
    !
    call Timers_stop("rad_resetData")
    !
  end subroutine rad_resetLWContainers




END SUBROUTINE RayTrace_LW



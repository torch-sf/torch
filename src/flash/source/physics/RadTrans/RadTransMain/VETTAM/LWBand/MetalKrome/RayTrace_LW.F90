!===============================================================================
!
! Subroutine: RayTrace_LW
! Path: source/physics/RadTrans/RadTransMain/VETTAM/LWBand/MetalKrome/RayTrace_LW.F90
!
!===============================================================================
!
! Description: Ray-Trace step for computing the dissociation of H2, CO and the ionization of C via LW photons. 
! Author: Shyam Menon (2024)
! Email: shyam.menon@rutgers.edu / smenon@flatironsinstitute.org
!


#undef RT_CORRECT
#undef DEBUG_RT


#include "Flash.h"
#include "constants.h"
subroutine RayTrace_LW()
  use RadTrans_RayTrace_3DRT
  use rt_lwdata
  use rt_lwmodule
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

  !The number of optical depth/column density calculations done in this ray-trace
  !5 species: Dust taud, and H2, CO, C, H column densities for the shielding factors
  integer, parameter  :: NLW_SPEC = 5
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
  real :: ColDens, ColDens_local, ColDens_global, taud
  real :: intensity_global
  real :: dest(NDIM)

  integer, dimension(:), allocatable :: id_sorted, QSindex
  !Containers to store the column density of H2 and the dust optical depth, needed for other species
  !Containers would store these for every cell in the block, for every particle -- i.e. shape (nAllParticles, NXB, NYB, NZB)
  real, dimension(:,:,:,:,:), allocatable :: NH2_star, taud_star
  ! Particle stuff
  integer              :: nAllParticles
  !Flag for which species the shielding is computed; 0: H2, 1: H
  character(len=MAX_STRING_LENGTH) :: fsh_spec



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

  ! allocate containers for NH2 and taud
  allocate(NH2_star(nrOfLeafBlocks,NXB,NYB,NZB,nAllParticles), taud_star(nrOfLeafBlocks,NXB,NYB,NZB,nAllParticles))

  fsh_spec = "dust" !Indicate that dust optical depth to be calculated
  call raytrace_species()

  fsh_spec = "H2" !Indicate that shielding function for H2 is to be done
  call raytrace_species()  

  !Compute the LW-flux weighted average shielding factor for H2
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do

  if(useHIshield) then
    !Now do the same steps for H
    fsh_spec = "H" !Indicate that shielding function for H is to be done
    call raytrace_species()
    !Compute the LW-flux weighted average shielding factor for H
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call average_fshield(blk)
    end do
  endif

#ifdef C_SPEC
  !Now do the same steps for C
  fsh_spec = "C" !Indicate that shielding function for C is to be done
  call raytrace_species()
  !Compute the LW-flux weighted average shielding factor for C
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do
#endif

#ifdef CO_SPEC
  !Now do the same steps for CO
  fsh_spec = "CO" !Indicate that shielding function for CO is to be done
  call raytrace_species()
  !Compute the LW-flux weighted average shielding factor for CO
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do
#endif

  !Compute the LW dissociation rate, C ionization rate and CO dissociation
  ! if(lwdiss_type .eq. "raytrace") then
  !   fsh_spec = "flux" !This is done to ensure only LWFL_VAR is reset in next step
  !   call rad_resetLWContainers()
  !   !Now combine above shielding factors to get the LW flux
  !   call calculate_LWFlux()
  ! endif

  deallocate(id_sorted,QSindex)
  deallocate(NH2_star, taud_star)
    

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Raytracing for Lyman Werner Radiation completed", "[RT_SINK] ")

  call Timers_stop("rad_raytrace_sink")

  CONTAINS

  subroutine raytrace_species()
    !This subroutine is used to perform the ray-tracing for each species column density
    !Steps: i) reset containers where flux-weighted sums are stored, ii) Set OPAC_VAR to the appropriate quantity for the species (noDens for column calculations, or kappa for taud)
    !       iii) Compute face contribution, communicate face values, and combine the SC solve using face values with LC solve to get the total column density

    !Reset the containers
    call rad_resetLWContainers()

    !Set OPAC_VAR
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


  end subroutine raytrace_species

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
    !The approach is to accumulate the total (dust-attenuated) flux in variable LWFL_VAR
    !The (dust-attenuated) flux-weighted sum over all sources for the shielding factors are stored in FSHM, FSHA, IC, DICO
    !The combination of the two is used to compute the effective shielfing factor for each species (this is done later)
    real :: luminosity, localColDens, rho, del(NDIM), radius_sink, flux, fshield
    real :: dColDens, dlcell, NH2
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

                if(fsh_spec == "dust") then
                  !For dust, no shielding factor calculation. Just store the optical depth for later use. 
                  taud_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle)) = ColDens
                  !return because no shielding calculations required here
                  cycle
                else
                  !Use the stored value to account for attenuation due to dust
                  taud = taud_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle))
                endif

                if(fsh_spec == "H2") then
                  !For dust, no shielding factor calculation. Just store the optical depth for later use. 
                  NH2_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle)) = ColDens
                else if(fsh_spec == "C" .or. fsh_spec == "CO") then
                  !Use the stored value to account for attenuation due to dust
                  NH2 = NH2_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle))
                endif

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


                !Attenuate flux due to dust absorption
                flux = flux * exp(-taud)

                !Get the shielding factor for the obtained column density
                if(fsh_spec == "H2") then
                  fshield = get_fshield_H2(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(FSHM_VAR,i,j,k) = solnData(FSHM_VAR,i,j,k) + flux * fshield
                else if(fsh_spec == "H") then
                  fshield = get_fshield_H(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(FSHA_VAR,i,j,k) = solnData(FSHA_VAR,i,j,k) + flux * fshield

                else if(fsh_spec == "C") then
                  !Obtain the shielding due to H2 and C for C
                  fshield = get_fshield_C(ColDens, NH2, bfive)
                  solnData(IC_VAR,i,j,k) = solnData(IC_VAR,i,j,k) + flux * fshield

                else if(fsh_spec == "CO") then
                  !Obtain the shielding due to H2 and CO for CO
                  fshield = get_fshield_CO(ColDens, NH2, bfive)
                  solnData(DICO_VAR,i,j,k) = solnData(DICO_VAR,i,j,k) + flux * fshield
                endif

                !Flux-weighted shielding: this is not the final value of the LW flux used
                !This counter is used to compute the average flux-weighted fshield over all species; will be replaced later
                solnData(LWFL_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) + flux 
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
    real :: nH2, nH, nC, nCO
    integer :: i,j,k
    real, dimension(NDIM) :: position, sink_pos

    call Grid_getBlkPtr(blk, solnData)
    do k=GRID_KLO_GC,GRID_KHI_GC
      do j=GRID_JLO_GC,GRID_JHI_GC
        do i=GRID_ILO_GC,GRID_IHI_GC
            
            if(fsh_spec == "H2") then
              !Obtain the number density of nH2
#ifdef H2_SPEC
              nH2 = solnData(DENS_VAR,i,j,k) * solnData(H2_SPEC,i,j,k)/h2A
#else
              call Driver_abortFlash("The Multispecies container H2_SPEC has to be present for LW dissociation.")
#endif
              !Now set the container OPAC_VAR which is the quantity integrated in the ray-tracer to the number density;
              !this ensures the column density is returned
              solnData(OPAC_VAR,i,j,k) = nH2
            else if(fsh_spec == "H") then
              !Obtain the number density of nH
#ifdef IHA_SPEC
              nH = solnData(DENS_VAR,i,j,k) * solnData(IHA_SPEC,i,j,k)/hA
#else
              call Driver_abortFlash("The Multispecies container IHA_SPEC has to be present for LW dissociation.")
#endif
              !Now set the container OPAC_VAR which is the quantity integrated in the ray-tracer to the number density;
              !this ensures the column density is returned
              solnData(OPAC_VAR,i,j,k) = nH
            else if(fsh_spec == "C") then
            !Obtain the number density of C
#ifdef C_SPEC
              nC = solnData(DENS_VAR,i,j,k) * solnData(C_SPEC,i,j,k)/cA
            
#else
              call Driver_abortFlash("The Multispecies container C_SPEC has to be present for C ionization.")
#endif
              solnData(OPAC_VAR,i,j,k) = nC
            else if(fsh_spec == "CO") then
            !Obtain the number density of CO
#ifdef CO_SPEC
              nCO = solnData(DENS_VAR,i,j,k) * solnData(CO_SPEC,i,j,k)/coA
#else
              call Driver_abortFlash("The Multispecies container CO_SPEC has to be present for CO dissociation.")
#endif
              solnData(OPAC_VAR,i,j,k) = nCO

            else if(fsh_spec == "dust") then
              !Obtain the dust opacity per unit length; this has already been set for the LW band in Cool.F90
              solnData(OPAC_VAR,i,j,k) = solnData(TAUP_VAR,i,j,k)
            endif

            
          

        end do 
      end do
    end do
    call Grid_releaseBlkPtr(blk, solnData)

  END SUBROUTINE SetOPAC_LW

  SUBROUTINE average_fshield(blk)
    !This subroutine computes the dust-attenuate flux-weighted shielding factor for H2, H, C and CO
    IMPLICIT NONE
    integer, intent(in) :: blk
    integer :: i,j,k
    real :: totalflux, fshield

    call Grid_getBlkPtr(blk, solnData)
    do k=GRID_KLO_GC,GRID_KHI_GC
        do j=GRID_JLO_GC,GRID_JHI_GC
            do i=GRID_ILO_GC,GRID_IHI_GC
                totalflux = solnData(LWFL_VAR,i,j,k)
                !The total cumulative flux from sinks is stored in FSHM_VAR or FSHA_VAR (depending on which is being updated now)
                if(fsh_spec == "H2") then
                  fshield = solnData(FSHM_VAR,i,j,k)
                else if(fsh_spec == "H") then
                  fshield = solnData(FSHA_VAR,i,j,k)
                else if(fsh_spec == "C") then
                  fshield = solnData(IC_VAR,i,j,k)
                else if(fsh_spec == "CO") then
                  fshield = solnData(DICO_VAR,i,j,k)
                endif
                !Now replace this var with the appropriate data, i.e. the avg shielding factor in FSHL_VAR
                if(totalflux .gt. 0) then
                    fshield = fshield/totalflux
                else
                    fshield = 0.0
                endif
                
                !Now store the shielding quantities in their containers
                if(fsh_spec == "H2") then
                  solnData(FSHM_VAR,i,j,k) = fshield
                else if(fsh_spec == "H") then
                  solnData(FSHA_VAR,i,j,k) = fshield
                else if(fsh_spec == "C") then
                  solnData(IC_VAR,i,j,k) = fshield
                else if(fsh_spec == "CO") then
                  solnData(DICO_VAR,i,j,k) = fshield
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
      if(fsh_spec == "H2") solnData(FSHM_VAR,:,:,:) = 0.0
      if(fsh_spec == "H") solnData(FSHA_VAR,:,:,:) = 0.0
      if(fsh_spec == "C") solnData(IC_VAR,:,:,:) = 0.0
      if(fsh_spec == "CO") solnData(DICO_VAR,:,:,:) = 0.0
      solnData(LWFL_VAR,:,:,:) = 0.0
      call Grid_releaseBlkPtr(blk, solnData)

    enddo ! b
    !
    call Timers_stop("rad_resetData")
    !
  end subroutine rad_resetLWContainers




END SUBROUTINE RayTrace_LW


